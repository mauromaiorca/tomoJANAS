[Repository home](../README.md) · [Documentation index](index.md) · [Accessory utilities](accessory_utils.md) · [CLI reference](reference/cli.md)

---

# sigma_estimate

`sigma_estimate` derives a suitable value of the Gaussian scale parameter σ used by the Structural Cross-correlation Index (SCI) directly from a pair of reference half-maps. The estimated σ can then be passed to JANAS scoring and session-manager commands via `--sigma`.

## Usage

```bash
janas_utils sigma_estimate halfmap1.mrc halfmap2.mrc [--mask mask.mrc]
```

The command prints the recommended σ in pixel units. If a 3D mask is supplied, the FSC is computed within the masked region so that σ reflects the local resolution of the region of interest.

## Background

In JANAS, SCI is built from first- and second-order spatial derivatives of amplitude-equalised images. Derivatives are computed in scale space, by convolution with a Gaussian of standard deviation σ. The choice of σ controls a trade-off:

- **Larger σ** increases regularisation of the derivative estimates and suppresses high-frequency noise, at the cost of sensitivity to fine structure.
- **Smaller σ** increases locality and sensitivity to small-scale features, at the cost of more noise amplification.

σ is specified in **pixel units**, so the corresponding physical scale is

$$\sigma_{\mathrm{\AA}} = \sigma \cdot \mathrm{apix}$$

where `apix` is the pixel size in Å/pixel. If particle images are resampled (for example by binning), σ should be adjusted to preserve the same physical scale.

σ can be selected manually or estimated automatically. `sigma_estimate` implements the automatic procedure described below.

## Automated estimation procedure

The estimator is FSC-linked: it derives σ from the Fourier Shell Correlation between the two reference half-maps, computed at the same pixel size used for scoring.

### Step 1 — Target frequency

The FSC at the gold-standard threshold 0.143 yields the spatial frequency *f*<sub>0.143</sub> in cycles/Å. JANAS caps the target frequency just below Nyquist to avoid edge artefacts:

$$f_{\mathrm{tgt}} = \min\!\bigl(f_{0.143},\; 0.95\, f_{\mathrm{Nyq}}\bigr)$$

where *f*<sub>Nyq</sub> = 0.5 / apix.

### Step 2 — Express target on the pixel grid

The target frequency is converted to cycles per pixel:

$$k_{\mathrm{tgt}} = f_{\mathrm{tgt}} \cdot \mathrm{apix} \quad \text{(cycles/pixel)}$$

### Step 3 — Raw σ from the Gaussian attenuation criterion

A raw estimate σ<sub>raw</sub> is chosen so that the Gaussian factor exp(−½ (2π σ *k*<sub>tgt</sub>)²) reaches a user-defined attenuation γ at *k*<sub>tgt</sub>, scaled by a user-controlled multiplier σ<sub>scale</sub>:

$$\sigma_{\mathrm{raw}} = \sigma_{\mathrm{scale}} \cdot \sqrt{\frac{2 \ln(1/\gamma)}{(2\pi)^2\, k_{\mathrm{tgt}}^2}}$$

Defaults: **γ = 0.5**, **σ<sub>scale</sub> = 2.0**.

### Step 4 — Frequency-dependent regularisation

When *f*<sub>0.143</sub> is low (and therefore *k*<sub>tgt</sub> is small), σ<sub>raw</sub> can become large and sensitive to small fluctuations in *f*<sub>0.143</sub>. JANAS applies a frequency-dependent regularisation towards a prior scale σ<sub>0</sub>, yielding the final estimate:

$$\sigma = \sigma_0 \left(\frac{\sigma_{\mathrm{raw}}}{\sigma_0}\right)^{k_{\mathrm{tgt}} / (k_{\mathrm{tgt}} + k_0)}$$

Defaults: **σ<sub>0</sub> = 1.0 pixel**, **k<sub>0</sub> = 0.20 cycles/pixel**.

*k*<sub>0</sub> sets the transition frequency at which the estimate moves from being dominated by the prior σ<sub>0</sub> to being dominated by σ<sub>raw</sub>.

### Multiple half-map pairs and masks

When more than one half-map pair is provided, σ is computed per pair and averaged. If a 3D mask is supplied, the FSC is computed within the masked region (after thresholding at 0.5, with optional softening by a user-defined rim width) so that σ reflects the local region of interest rather than the whole map.

## Choosing σ manually

For manual selection, σ is treated as a scale parameter controlling the trade-off described above. A practical approach is to run JANAS with a small set of σ values and compare the downstream reconstructions under otherwise identical processing conditions.

## Reference

The full derivation and rationale are described in the JANAS manuscript (Methods section: "Choice of the scale parameter σ"). See [Citation](citation.md).

---

[Back to Accessory utilities](accessory_utils.md) · [Back to documentation index](index.md)
