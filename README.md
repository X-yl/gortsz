# gortsz: gorilla time series compression

This is a `#![no_std]` and no `alloc` crate implementing Facebook's [gorilla](https://www.vldb.org/pvldb/vol8/p1816-teller.pdf)
compression algorithm for time-series data.

Basic usage:

```rs
use gortsz::*;
use gortsz::stats::*;

fn main() {
    let data = [
        (10u32, [1.0, 2.0, 3.3]),
        (11u32, [1.1, 2.2, 2.2]),
        (14u32, [1.2, 2.3, 3.8]),
        (16u32, [1.3, 2.9, 4.9]),
        (18u32, [1.4, 1.0, 3.2]),
    ];

    let buf = [0u8; 256];
    let compressed = match compress::<3, 3, WhitepaperOptions>(
        timeseries.iter(),
        buf.view_bits_mut(),
    ) {
        Ok(compressed) => compressed,
        Err(CompressError {
            valid_bits, ..
        }) => {
            valid_bits // Buffer too small, only partly compressed
        }
    };

    for result in decompressor {
        match result {
            Ok((time, data)) => {
                println!("Decompressed time: {}, data: {:?}", time, data);
                println!("Original time: {}, data: {:?}", reference.0, reference.1);
                assert!(time == reference.0);
                assert!(data == reference.1);
            }
            Err(e) => {
                eprintln!("Decompression error: {}", e);
            }
        }
    }
}
```

The `WhitepaperOptions` struct implements the trait `CompressionOptions` which gives
the options described in the whitepaper. Implement this trait to use different options.