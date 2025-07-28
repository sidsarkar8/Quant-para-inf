# Inference for Quantile-Parametrized Distributions

This repository contains code and methods for developing inference procedures tailored to **quantile-based distributions**, which are flexible yet analytically challenging.

## Motivation

Quantile-parametrized distributions (e.g., Tukey's Lambda, Generalized Lambda Distribution) are widely used due to their ability to capture a variety of shapes and tail behaviors. However, standard inference techniques encounter major roadblocks:

- ❌ **No closed-form density**: Prevents likelihood-based methods from being directly applied.
- ❌ **Irregular asymptotics**: Maximum Likelihood Estimation (MLE) may exhibit **non-$\sqrt{n}$ convergence rates** and **non-normal limiting distributions**, especially in certain regions of the parameter space.
- ❌ **Unreliable resampling**: Bootstrap and related methods may fail to approximate the true sampling distribution.

## Our Contribution

We develop a **new inference framework** specifically designed for quantile-parametrized models. The key features are:

- ✅ **Quantile-function based**: Inference is performed using the quantile representation directly, avoiding density estimation.
- ✅ **Assumption-lean**: Requires minimal structural assumptions, allowing greater robustness.
- ✅ **Principled guarantees**: The method retains statistical validity even when classical assumptions fail.

## Applications

- Inference for location, scale, or shape parameters in quantile-based families.
- Use in simulation-based models where only quantile functions are available.
- Robust alternatives in applied settings where tail behavior is complex or unknown.


