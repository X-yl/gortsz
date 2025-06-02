#![no_std]
#![feature(generic_arg_infer)]

use core::{cell::Cell, error::Error, fmt::Display, ops::Not};

pub use bitvec::prelude::*;
use stats::{CompressionOptions, CompressionStatsExt, NoStats};

pub mod stats;

// sign contract for i32.
// NB: this returns a BitArray which will be u8 aligned. Callers must
//     slice the first $n bits to get the correct result.
macro_rules! sign_contract_32 {
    ($value:expr, $n:expr) => {{
        let mut result = bitarr!(u8, Lsb0; 0; 32);
        let value = $value as i32;
        value
            .to_le_bytes()
            .view_bits::<Lsb0>()
            .into_iter()
            .enumerate()
            .filter_map(|(i, b)| if i < $n - 1 || i == 31 { Some(b) } else { None })
            .enumerate()
            .for_each(|(i, b)| result.set(i, *b));
        result
    }};
}

pub struct CompressError<'a> {
    pub valid_bits: &'a BitSlice<u8>,
    pub entries_processed: usize,
}

pub fn compress<'a, 'b, const N: usize, const BC: usize, OPTS: CompressionOptions<BC>>(
    series: impl IntoIterator<Item = &'a (u32, [f32; N])>,
    buf: &'b mut BitSlice<u8>,
) -> Result<&'b BitSlice<u8>, CompressError<'b>> {
    compress_with_stats(series, buf, &mut NoStats::<BC, OPTS>::new())
}
/// Facebook's Gorilla compression algorithm for time series data.
/// The timestamps are delta-to-delta encoded and the data is XOR encoded.
pub fn compress_with_stats<
    'a,
    'b,
    const N: usize,
    const BC: usize,
    OPTS: CompressionOptions<BC>,
>(
    series: impl IntoIterator<Item = &'a (u32, [f32; N])>,
    buf: &'b mut BitSlice<u8>,
    stats: &mut impl CompressionStatsExt<BC, OPTS>,
) -> Result<&'b BitSlice<u8>, CompressError<'b>> {
    let mut series = series.into_iter();
    let mut index = 0;
    let last_valid_index = Cell::new(0);

    macro_rules! write_bits_ {
        ($b:expr, $idx:expr) => {
            if buf.len() < index + $b.len() {
                return Err(CompressError {
                    valid_bits: &buf[..last_valid_index.get()],
                    entries_processed: $idx,
                });
            }
            buf[index..(index + $b.len())].copy_from_bitslice($b);
            index += $b.len();
        };
    }
    macro_rules! write_bit_ {
        ($b:expr, $idx:expr) => {
            if buf.len() < index + 1 {
                return Err(CompressError {
                    valid_bits: &buf[..last_valid_index.get()],
                    entries_processed: $idx,
                });
            }
            buf.set(index, $b);
            index += 1;
        };
    }

    let (mut previous_time, mut previous_data) =
        *series.next().expect("Time series must not be empty");
    let mut previous_delta = 0i32;
    let mut previous_xors = [0u32; N];

    write_bits_!(previous_time.to_le_bytes().view_bits::<Lsb0>(), 0);
    for bytes in previous_data.iter().map(|it| it.to_le_bytes()) {
        write_bits_!(BitSlice::<u8, Lsb0>::from_slice(&bytes), 0);
    }

    for (row_index, &(time, ref data)) in series.enumerate() {
        macro_rules! write_bit {
            ($b:expr) => {
                write_bit_!($b, row_index);
            };
        }
        macro_rules! write_bits {
            ($b:expr) => {
                write_bits_!($b, row_index);
            };
        }

        let delta = (time as i64 - previous_time as i64) as i32;
        let delta_delta = delta - previous_delta;

        if delta_delta == 0 {
            write_bits!(bits![u8, Lsb0; 0u8]);
            stats.increment_repeated_count();
        } else {
            let mut fit_in_bin = false;
            for (i, bin) in OPTS::DELTA_DELTA_BINS.iter().enumerate() {
                if delta_delta >= -(1 << (bin - 1)) && delta_delta < (1 << (bin - 1)) {
                    // Write i + 1 one bits then a zero bit
                    for _ in 0..(i + 1) {
                        write_bits!(bits![u8, Lsb0; 1]);
                    }
                    write_bits!(bits![u8, Lsb0; 0]);

                    // Write the delta delta
                    let bits = &sign_contract_32!(delta_delta, *bin as usize)[..*bin as usize];
                    write_bits!(&bits);
                    fit_in_bin = true;
                    stats.increment_bin(i);
                    break;
                }
            }
            if !fit_in_bin {
                // Write n+1 one bits and no zero bit
                for _ in 0..(OPTS::DELTA_DELTA_BINS.len() + 1) {
                    write_bits!(bits![u8, Lsb0; 1]);
                }

                // Write the delta delta
                let bits = sign_contract_32!(delta_delta, 32);
                stats.increment_overflow_count();
                write_bits!(&bits);
            }
        }

        previous_delta = delta;
        previous_time = time;

        for ((&d, previous_d), previous_xor) in data
            .iter()
            .zip(previous_data.iter_mut())
            .zip(previous_xors.iter_mut())
        {
            let xor = d.to_bits() ^ previous_d.to_bits();

            if xor == 0 {
                // No change, write a zero bit
                write_bits!(bits![u8, Lsb0; 0]);
                *previous_xor = 0;
                continue;
            } else {
                write_bits!(bits![u8, Lsb0; 1]);
            }

            let leading_zeros = xor.leading_zeros();
            let trailing_zeros = xor.trailing_zeros();
            let prev_leading_zeros = previous_xor.leading_zeros();
            let prev_trailing_zeros = previous_xor.trailing_zeros();

            // NB: if previous_xor is zero, we will read 32 for prev_leading_zeros and prev_trailing_zeros
            //     however, this is fine as in the two cases:
            //      - xor is also zero, in which case we already handled it above
            //      - xor is not zero, in which case leading_zeros and trailing_zeros will be less than 32
            //        so we will enter the second branch below.

            if leading_zeros >= prev_leading_zeros && trailing_zeros >= prev_trailing_zeros {
                // The new data is fully contained within the previous data's leading and trailing zeros
                write_bits!(bits![u8, Lsb0; 0]);

                let n_bits_to_write = 32 - prev_leading_zeros - prev_trailing_zeros;
                let compressed = xor
                    .view_bits::<Lsb0>()
                    .iter()
                    .skip(prev_trailing_zeros as usize)
                    .take(n_bits_to_write as usize);

                for bit in compressed {
                    let bit = *bit;
                    write_bit!(bit);
                }
            } else {
                let leading_zeros = leading_zeros.min(0b1111);

                write_bits!(bits![u8, Lsb0; 1]);

                // this is kind of a weird unexplained decision in the whitepaper,
                // why use 5 bits for a f64 (adapted to 4 bits for f32 here)?
                let leading_zero_bits = leading_zeros.to_le_bytes();
                write_bits!(&leading_zero_bits.view_bits::<Lsb0>()[..4]);

                let n_meaningful_bits = 32 - leading_zeros - trailing_zeros;
                write_bits!(&(n_meaningful_bits - 1).to_le_bytes().view_bits::<Lsb0>()[..5]);

                for b in xor
                    .view_bits::<Lsb0>()
                    .iter()
                    .skip(trailing_zeros as usize)
                    .take(n_meaningful_bits as usize)
                {
                    write_bit!(*b);
                }
            }

            *previous_d = d;
            *previous_xor = xor;
        }

        last_valid_index.set(index);
    }

    Ok(buf.get(0..index).expect("checked bounds already"))
}

pub struct Decompressor<'a, const N: usize, const BC: usize, OPTS: CompressionOptions<BC>> {
    buf: &'a BitSlice<u8>,
    index: usize,
    previous_time: u32,
    previous_delta: i32,
    previous_data: [f32; N],
    previous_xors: [u32; N],
    failed: bool,
    _options: core::marker::PhantomData<OPTS>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompressError {
    MissingHeader,
    CorruptedTimestamp { index: usize },
    CorruptedData { index: usize },
}

impl Display for DecompressError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DecompressError::MissingHeader => write!(f, "Missing header in the compressed data"),
            DecompressError::CorruptedTimestamp { index } => {
                write!(f, "Corrupted timestamp at bit {index}")
            }
            DecompressError::CorruptedData { index } => {
                write!(f, "Corrupted data at bit {index}")
            }
        }
    }
}

impl Error for DecompressError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        None
    }

    fn description(&self) -> &str {
        "description() is deprecated; use Display"
    }

    fn cause(&self) -> Option<&dyn Error> {
        self.source()
    }
}

impl<'a, const N: usize, const BC: usize, OPTS: CompressionOptions<BC>>
    Decompressor<'a, N, BC, OPTS>
{
    pub fn new(buf: &'a BitSlice<u8>) -> Self {
        Self {
            buf,
            index: 0,
            previous_time: 0,
            previous_data: [0.0; N],
            previous_delta: 0,
            previous_xors: [0; N],
            failed: false,
            _options: core::marker::PhantomData,
        }
    }
}

impl<'a, const SAMPLES_PER_ROW: usize, const BC: usize, OPTS: CompressionOptions<BC>> Iterator
    for Decompressor<'a, SAMPLES_PER_ROW, BC, OPTS>
{
    type Item = Result<(u32, [f32; SAMPLES_PER_ROW]), DecompressError>;

    fn next(&mut self) -> Option<Self::Item> {
        let consume = |index: &mut usize, n: usize| {
            if self.buf.len() - *index < n {
                return None;
            }
            let bits = self.buf.get(*index..*index + n)?;
            *index += n;
            Some(bits)
        };
        if self.failed || self.index >= self.buf.len() {
            return None; // No more data to decompress
        }

        if self.index == 0 {
            if self.buf.len() < 32 + (SAMPLES_PER_ROW * 32) {
                self.failed = true;
                return Some(Err(DecompressError::MissingHeader));
            }

            let first_time: u32 = consume(&mut self.index, 32)
                .expect("checked bounds already")
                .load_le();
            let first_data: [f32; SAMPLES_PER_ROW] = self.buf[32..(32 + (SAMPLES_PER_ROW * 32))]
                .chunks_exact(32)
                .map(|chunk| f32::from_bits(chunk.load_le()))
                .collect::<heapless::Vec<f32, SAMPLES_PER_ROW>>()
                .into_array()
                .expect("Impossible for N f32's to not fit into an array of N f32's");

            self.index = 32 + (SAMPLES_PER_ROW * 32);
            self.previous_time = first_time;
            self.previous_data = first_data;

            return Some(Ok((self.previous_time, self.previous_data)));
        }

        let timestamp_ctl_bits = self
            .buf
            .get(self.index..)
            .expect("Checked bounds already")
            .leading_ones();

        self.index += core::cmp::min(timestamp_ctl_bits + 1, OPTS::DELTA_DELTA_BINS.len() + 1);

        if timestamp_ctl_bits == 0 {
            // No change in delta
            if let Some(t) = self.previous_time.checked_add_signed(self.previous_delta) {
                self.previous_time = t;
            } else {
                self.failed = true;
                return Some(Err(DecompressError::CorruptedTimestamp {
                    index: self.index,
                }));
            }
        } else {
            let n_dd_bits = *OPTS::DELTA_DELTA_BINS
                .get(timestamp_ctl_bits - 1)
                .unwrap_or(&32) as usize;

            let delta_delta: i32 = if let Some(dd) = consume(&mut self.index, n_dd_bits) {
                dd.load_le()
            } else {
                self.failed = true;
                return Some(Err(DecompressError::CorruptedTimestamp {
                    index: self.index,
                }));
            };

            let delta = self.previous_delta + delta_delta;
            self.previous_time = self.previous_time.wrapping_add(delta as u32);
            self.previous_delta = delta;
        };

        for (previous_d, previous_xor) in self
            .previous_data
            .iter_mut()
            .zip(self.previous_xors.iter_mut())
        {
            let new_data_bit = if let Some(x) = consume(&mut self.index, 1) {
                x[0]
            } else {
                self.failed = true;
                return Some(Err(DecompressError::CorruptedTimestamp {
                    index: self.index,
                }));
            };

            if new_data_bit.not() {
                // No change, just copy the previous value
                *previous_xor = 0;
                continue;
            } else {
                let xor_ctl_bit = if let Some(x) = consume(&mut self.index, 1) {
                    x[0]
                } else {
                    self.failed = true;
                    return Some(Err(DecompressError::CorruptedData { index: self.index }));
                };

                if xor_ctl_bit.not() {
                    // The new data is fully contained within the previous data's leading and trailing zeros
                    let n_bits_to_read =
                        (32 - previous_xor.leading_zeros() - previous_xor.trailing_zeros())
                            as usize;
                    let compressed = if let Some(x) = consume(&mut self.index, n_bits_to_read) {
                        x
                    } else {
                        self.failed = true;
                        return Some(Err(DecompressError::CorruptedData { index: self.index }));
                    };

                    let xor = compressed.load_le::<u32>() << previous_xor.trailing_zeros();
                    *previous_d = f32::from_bits(xor ^ previous_d.to_bits());
                    *previous_xor = xor;
                } else {
                    let leading_zero_bits = if let Some(x) = consume(&mut self.index, 4) {
                        x.load_le::<u32>()
                    } else {
                        self.failed = true;
                        return Some(Err(DecompressError::CorruptedData { index: self.index }));
                    };

                    let meaningful_bits_len = if let Some(x) = consume(&mut self.index, 5) {
                        // + 1 to convert to 0-based indexing
                        x.load_le::<u32>() + 1
                    } else {
                        self.failed = true;
                        return Some(Err(DecompressError::CorruptedData { index: self.index }));
                    };

                    if leading_zero_bits + meaningful_bits_len > 32 {
                        self.failed = true;
                        return Some(Err(DecompressError::CorruptedData { index: self.index }));
                    }

                    let trailing_zeros = 32 - leading_zero_bits - meaningful_bits_len;

                    let meaningful_bits =
                        if let Some(x) = consume(&mut self.index, meaningful_bits_len as usize) {
                            x
                        } else {
                            self.failed = true;
                            return Some(Err(DecompressError::CorruptedData { index: self.index }));
                        };

                    let xor = meaningful_bits.load_le::<u32>() << trailing_zeros;
                    *previous_d = f32::from_bits(xor ^ previous_d.to_bits());
                    *previous_xor = xor;
                }
            }
        }

        return Some(Ok((self.previous_time, self.previous_data)));
    }
}

#[cfg(test)]
extern crate std;

#[cfg(test)]
mod tests {
    use crate::stats::{CompressionStats, whitepaper::WhitepaperOptions};
    use std::{eprintln, println, vec::Vec};

    use super::*;

    #[test]
    fn compression_test() {
        const TIMESERIES: &'static str = include_str!("../samples.csv");
        const COLUMNS: usize = const_str::split!(const_str::split!(TIMESERIES, '\n')[0], ',').len();

        let timeseries: Vec<_> = TIMESERIES
            .lines()
            .map(|line| {
                let mut spl = line.split(',');
                let time = spl.next().unwrap().trim().parse::<f32>().unwrap() as u32;
                let samples: [f32; COLUMNS - 1] = spl
                    .map(|col| col.trim().parse::<f32>().unwrap())
                    .collect::<heapless::Vec<_, { COLUMNS - 1 }>>()
                    .into_array()
                    .expect("checked size already");
                (time, samples)
            })
            .collect();

        let mut buf = [0u8; 250];
        let mut entries_processed = timeseries.len();
        let compressed = match compress::<{ COLUMNS - 1 }, 3, WhitepaperOptions>(
            timeseries.iter(),
            buf.view_bits_mut(),
        ) {
            Ok(compressed) => compressed,
            Err(CompressError {
                valid_bits,
                entries_processed: entries,
            }) => {
                eprintln!(
                    "warn: Buffer too small, only compressed {entries} entries out of {}",
                    timeseries.len()
                );
                entries_processed = entries;
                valid_bits
            }
        };

        let compressed_size_bits = compressed.len();
        let data_size_bytes = entries_processed * (4 + 4 * timeseries[0].1.len());

        let decompressor = Decompressor::<_, _, WhitepaperOptions>::new(compressed);
        for (result, reference) in decompressor.zip(timeseries.iter()) {
            match result {
                Ok((time, data)) => {
                    assert!(time == reference.0);
                    assert!(data == reference.1);
                }
                Err(e) => {
                    panic!("Decompression failed: {e}");
                }
            }
        }

        println!(
            "Compressed size: {} bytes, originally {} bytes. Compression ratio: {:.2}%",
            compressed_size_bits / 8,
            data_size_bytes,
            compressed_size_bits as f64 / 8. * 100.0 / data_size_bytes as f64
        );
    }
}
