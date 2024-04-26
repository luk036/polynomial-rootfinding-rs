# 🫚 ginger-rs

Polynomial root-finding algorithms (parallelizable) in Rust

The Aberth-Ehrlich (AE) method introduced by Ehrlich (1967); Aberth (1973) combines the elements of Newton’s method with an implicit deflation strategy, which allows for the computation of all roots of a polynomial simultaneously and converges cubically. This method is considered an improvement over the Durand-Kerner method, another simultaneous root solver method which converges quadratically and is 10-100 times slower (see, e.g., Ghidouche et al. 2017). The facts that AE is extremely fast for various degrees of polynomials, its ability to find all the roots at once (unlike other iterative root-solving methods such as Laguerre’s) and its root polishing procedure, which is inside the main iteration loop and can be controlled for a chosen accuracy.

Bairstow's method is an iterative method that is used to find complex roots of a polynomial. This method is based on synthetic division and can be used to find all roots of a polynomial.

Parallel Bairstow's method refers to a parallel algorithm based on Bairstow's method. This algorithm uses parallel computation to accelerate the process of finding complex roots of polynomials.
[![Crates.io](https://img.shields.io/crates/v/ginger-rs.svg)](https://crates.io/crates/ginger-rs)
[![Docs.rs](https://docs.rs/ginger-rs/badge.svg)](https://docs.rs/ginger-rs)
[![CI](https://github.com/luk036/ginger-rs/workflows/CI/badge.svg)](https://github.com/luk036/ginger-rs/actions)
[![codecov](https://codecov.io/gh/luk036/ginger-rs/branch/main/graph/badge.svg?token=1qz6WD6Rs5)](https://codecov.io/gh/luk036/ginger-rs)

## 🛠️ Installation

### 📦 Cargo

- Install the rust toolchain in order to have cargo installed by following
  [this](https://www.rust-lang.org/tools/install) guide.
- run `cargo install ginger-rs`

## 👀 See also

- [ginger](https://luk036.github.io/ginger)
- [ginger-cpp](https://luk036.github.io/ginger-cpp)

## 📜 License

Licensed under either of

- Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license
  ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## 🤝 Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

See [CONTRIBUTING.md](CONTRIBUTING.md).
