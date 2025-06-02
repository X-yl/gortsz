/// A trait defining compression options.
pub trait CompressionOptions<const N: usize> {
    /// In order to compress the delta-delta values, we assign bins based on the
    /// number of bits required to represent the value.
    /// This should be an **ascending** array of number of bits that should fit in each bin.
    /// The whitepaper suggests [7, 9, 12].
    const DELTA_DELTA_BINS: [u8; N];
}

pub mod whitepaper {
    use super::CompressionOptions;

    #[derive(Debug, Clone, Copy)]
    pub struct WhitepaperOptions;
    impl CompressionOptions<3> for WhitepaperOptions {
        const DELTA_DELTA_BINS: [u8; 3] = [7, 9, 12];
    }
}

#[derive(Debug, Clone, Copy)]
pub enum XorCtl {
    Z,
    OZ,
    OO,
}

pub trait CompressionStatsExt<const N: usize, OPTS: CompressionOptions<N>> {
    fn increment_bin(&mut self, bin: usize);
    fn increment_repeated_count(&mut self);
    fn increment_overflow_count(&mut self);

    fn increment_data_xor_ctl(&mut self, xor_ctl: XorCtl);
}

#[derive(Debug, Clone, Copy)]
pub struct CompressionStats<const N: usize, OPTS: CompressionOptions<N>> {
    pub timestamp_delta_bin_distribution: [u64; N],
    pub repeated_time_count: u64,
    pub overflow_time_count: u64,

    pub xor_ctl0_count: u64,
    pub xor_ctl10_count: u64,
    pub xor_ctl11_count: u64,

    _options: core::marker::PhantomData<OPTS>,
}

pub struct NoStats<const N: usize, OPTS: CompressionOptions<N>> {
    _options: core::marker::PhantomData<OPTS>,
}

impl<const N: usize, OPTS: CompressionOptions<N>> NoStats<N, OPTS> {
    pub fn new() -> Self {
        Self {
            _options: core::marker::PhantomData,
        }
    }
}
impl<const N: usize, OPTS: CompressionOptions<N>> CompressionStatsExt<N, OPTS>
    for NoStats<N, OPTS>
{
    fn increment_bin(&mut self, _bin: usize) {}

    fn increment_repeated_count(&mut self) {}

    fn increment_overflow_count(&mut self) {}

    fn increment_data_xor_ctl(&mut self, _xor_ctl: XorCtl) {}
}

impl<const N: usize, OPTS: CompressionOptions<N>> CompressionStats<N, OPTS> {
    pub fn new() -> Self {
        Self {
            timestamp_delta_bin_distribution: [0; N],
            repeated_time_count: 0,
            overflow_time_count: 0,

            xor_ctl0_count: 0,
            xor_ctl10_count: 0,
            xor_ctl11_count: 0,

            _options: core::marker::PhantomData,
        }
    }
}

impl<const N: usize, OPTS: CompressionOptions<N>> CompressionStatsExt<N, OPTS>
    for CompressionStats<N, OPTS>
{
    fn increment_repeated_count(&mut self) {
        self.repeated_time_count += 1;
    }

    fn increment_overflow_count(&mut self) {
        self.overflow_time_count += 1;
    }

    fn increment_bin(&mut self, bin: usize) {
        self.timestamp_delta_bin_distribution[bin] =
            self.timestamp_delta_bin_distribution[bin].saturating_add(1);
    }

    fn increment_data_xor_ctl(&mut self, xor_ctl: XorCtl) {
        match xor_ctl {
            XorCtl::Z => self.xor_ctl0_count += 1,
            XorCtl::OZ => self.xor_ctl10_count += 1,
            XorCtl::OO => self.xor_ctl11_count += 1,
        }
    }
}
