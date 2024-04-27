# 🫚 polynomial-rootfinding-rs

Polynomial root-finding algorithms (parallelizable) in Rust

[![Crates.io](https://img.shields.io/crates/v/polynomial-rootfinding-rs.svg)](https://crates.io/crates/polynomial-rootfinding-rs)
[![Docs.rs](https://docs.rs/polynomial-rootfinding-rs/badge.svg)](https://docs.rs/polynomial-rootfinding-rs)
[![CI](https://github.com/luk036/polynomial-rootfinding-rs/workflows/CI/badge.svg)](https://github.com/luk036/polynomial-rootfinding-rs/actions)
[![codecov](https://codecov.io/gh/luk036/polynomial-rootfinding-rs/graph/badge.svg?token=q5sX82bm7S)](https://codecov.io/gh/luk036/polynomial-rootfinding-rs)

The Aberth-Ehrlich (AE) method introduced by Ehrlich (1967); Aberth (1973) combines the elements of Newton’s method with an implicit deflation strategy, which allows for the computation of all roots of a polynomial simultaneously and converges cubically. This method is considered an improvement over the Durand-Kerner method, another simultaneous root solver method which converges quadratically and is 10-100 times slower (see, e.g., Ghidouche et al. 2017). The facts that AE is extremely fast for various degrees of polynomials, its ability to find all the roots at once (unlike other iterative root-solving methods such as Laguerre’s) and its root polishing procedure, which is inside the main iteration loop and can be controlled for a chosen accuracy.

Bairstow's method is an iterative method that is used to find complex roots of a polynomial. This method is based on synthetic division and can be used to find all roots of a polynomial.

Parallel Bairstow's method refers to a parallel algorithm based on Bairstow's method. This algorithm uses parallel computation to accelerate the process of finding complex roots of polynomials.

## 🛠️ Installation

### 📦 Cargo

- Install the rust toolchain in order to have cargo installed by following
  [this](https://www.rust-lang.org/tools/install) guide.
- run `cargo install polynomial-rootfinding-rs`

## 👀 See also

- [polynomial_rootfinding](https://luk036.github.io/polynomial_rootfinding)
- [polynomial_rootfinding-cpp](https://luk036.github.io/polynomial_rootfinding-cpp)

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
