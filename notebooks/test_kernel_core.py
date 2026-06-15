#!/usr/bin/env python3
"""Verify MakeConvolutionKernel works without target PSF and can save processed PSFs.

Uses F335M and F770W as test cases.

Also compares the processed source PSF output with  
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits

from jwst_kernels.kernel_core import MakeConvolutionKernel
from kernel_core_original import (
    MakeConvolutionKernel as MakeConvolutionKernelOriginal,
)
from jwst_kernels.make_psf import read_PSF, save_nircam_PSF, save_miri_PSF

HERE = Path(__file__).parent
PSF_DIR = HERE / "psf_cache"
OUT_DIR = HERE


def main():
    PSF_DIR.mkdir(exist_ok=True)

    print(f"PSF cache: {PSF_DIR}")
    print(f"Output dir: {OUT_DIR}\n")

    print("Generating NIRCam F335M PSF...")
    save_nircam_PSF(["F335M"], output_dir=str(PSF_DIR))

    print("Generating MIRI F770W PSF...")
    save_miri_PSF(["F770W"], output_dir=str(PSF_DIR))

    source_data, source_pix = read_PSF(
        {"camera": "NIRCam", "filter": "F335M"}, psf_dir=str(PSF_DIR)
    )
    target_data, target_pix = read_PSF(
        {"camera": "MIRI", "filter": "F770W"}, psf_dir=str(PSF_DIR)
    )

    print(f"\nSource PSF shape: {source_data.shape}, pixscale: {source_pix}")
    print(f"Target PSF shape: {target_data.shape}, pixscale: {target_pix}")

    print("\n--- Test 1: process_source_psf with default options ---")
    ck1 = MakeConvolutionKernel(
        source_psf=source_data,
        source_pixscale=source_pix,
        source_name='F335M',
        verbose=True,
    )
    ck1.process_source_psf()
    print(f"Processed source shape: {ck1.source_psf.shape}")
    print(f"Resolved common_pixscale: {ck1.common_pixscale}")
    print(f"Resolved grid_size_arcsec: {ck1.grid_size_arcsec}")
    assert ck1.target_psf is None, "target_psf should be None"
    assert np.isclose(ck1.common_pixscale, source_pix), (
        "common_pixscale should default to source_pixscale"
    )

    saved = ck1.save_processed_psf(str(OUT_DIR), which='source')
    hdu = fits.open(saved[0])
    h = hdu[0].header
    print(f"  Saved: {saved[0]}")
    print(f"  shape={hdu[0].data.shape}  PSFNAME={h['PSFNAME']}")
    hdu.close()
    print("PASS\n")

    # ---- Test 2: make_convolution_kernel fails without target ----
    print("--- Test 2: make_convolution_kernel raises without target ---")
    try:
        kc2 = MakeConvolutionKernel(
            source_psf=source_data,
            source_pixscale=source_pix,
            common_pixscale=source_pix,
        )
        kc2.make_convolution_kernel()
        print("FAIL: should have raised ValueError")
    except ValueError as e:
        print(f"  Caught expected error: {e}")
        print("PASS\n")

    # ---- Test 3: full kernel creation still works ----
    print("--- Test 3: full kernel creation with both PSFs ---")
    common_pix = source_pix
    grid_size = np.array([361 * common_pix, 361 * common_pix])

    kc3 = MakeConvolutionKernel(
        source_psf=source_data,
        source_pixscale=source_pix,
        source_name='F335M',
        target_psf=target_data,
        target_pixscale=target_pix,
        target_name='F770W',
        common_pixscale=common_pix,
        grid_size_arcsec=grid_size,
        verbose=True,
    )
    kc3.make_convolution_kernel()
    print(f"Kernel shape: {kc3.kernel.shape}")

    saved2 = kc3.save_processed_psf(str(OUT_DIR), which='both')
    for p in saved2:
        print(f"  Saved: {p}")
    print("PASS\n")

    # ---- Test 4: parity with original pipeline ----
    print("--- Test 4: parity check vs original kernel_core ---")

    kc_orig = MakeConvolutionKernelOriginal(
        source_psf=source_data,
        source_pixscale=source_pix,
        source_name='F335M',
        target_psf=target_data,
        target_pixscale=target_pix,
        target_name='F770W',
        common_pixscale=common_pix,
        grid_size_arcsec=grid_size,
        verbose=True,
    )
    kc_orig.make_convolution_kernel()
    original_processed = kc_orig.source_psf_processed

    kc_new = MakeConvolutionKernel(
        source_psf=source_data,
        source_pixscale=source_pix,
        source_name='F335M',
        common_pixscale=common_pix,
        grid_size_arcsec=grid_size,
        verbose=True,
    )
    kc_new.process_source_psf()
    new_processed = kc_new.source_psf

    print(f"  Original shape: {original_processed.shape}")
    print(f"  New shape:      {new_processed.shape}")

    if original_processed.shape != new_processed.shape:
        print("FAIL: shapes differ!")
    else:
        max_abs_diff = np.max(np.abs(original_processed - new_processed))
        rms_diff = np.sqrt(np.mean((original_processed - new_processed) ** 2))
        max_val = np.max(np.abs(original_processed))
        print(f"  Max |diff|:           {max_abs_diff:.3e}")
        print(f"  RMS diff:             {rms_diff:.3e}")
        print(f"  Max |diff| / max val: {max_abs_diff / max_val:.3e}")

        if np.allclose(original_processed, new_processed, atol=1e-10):
            print("PASS: outputs are identical (atol=1e-10)\n")
        else:
            print("WARN: small numerical differences present\n")

    # ---- Plot ----
    fig, axes = pl.subplots(2, 3, figsize=(15, 9))

    im0 = axes[0, 0].imshow(
        np.log10(ck1.source_psf / ck1.source_psf.max() + 1e-8),
        origin="lower", cmap="inferno", vmin=-5, vmax=0,
    )
    axes[0, 0].set_title(f"F335M source (native grid)\n{ck1.source_psf.shape}")
    pl.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(
        np.log10(kc3.source_psf / kc3.source_psf.max() + 1e-8),
        origin="lower", cmap="inferno", vmin=-5, vmax=0,
    )
    axes[0, 1].set_title(f"F335M source (new)\n{kc3.source_psf.shape}")
    pl.colorbar(im1, ax=axes[0, 1])

    im_orig = axes[0, 2].imshow(
        np.log10(original_processed / original_processed.max() + 1e-8),
        origin="lower", cmap="inferno", vmin=-5, vmax=0,
    )
    axes[0, 2].set_title(f"F335M source (original)\n{original_processed.shape}")
    pl.colorbar(im_orig, ax=axes[0, 2])

    im2 = axes[1, 0].imshow(
        np.log10(kc3.target_psf / kc3.target_psf.max() + 1e-8),
        origin="lower", cmap="inferno", vmin=-5, vmax=0,
    )
    axes[1, 0].set_title(f"F770W target\n{kc3.target_psf.shape}")
    pl.colorbar(im2, ax=axes[1, 0])

    im3 = axes[1, 1].imshow(
        kc3.kernel / kc3.kernel.max(),
        origin="lower", cmap="RdBu_r", vmin=-0.05, vmax=1,
    )
    axes[1, 1].set_title(f"F335M->F770W kernel\n{kc3.kernel.shape}")
    pl.colorbar(im3, ax=axes[1, 1])

    if original_processed.shape == new_processed.shape:
        diff = new_processed - original_processed
        vlim = max(np.abs(diff.min()), np.abs(diff.max()))
        if vlim == 0:
            vlim = 1e-15
        im_diff = axes[1, 2].imshow(
            diff, origin="lower", cmap="RdBu_r", vmin=-vlim, vmax=vlim,
        )
        axes[1, 2].set_title("new - original")
        pl.colorbar(im_diff, ax=axes[1, 2])
    else:
        axes[1, 2].text(0.5, 0.5, "shape mismatch", ha="center", va="center",
                        transform=axes[1, 2].transAxes)
        axes[1, 2].set_title("new - original")

    fig.suptitle("kernel_core.py testing", fontsize=13)
    fig.subplots_adjust(wspace=0.4, hspace=0.35)

    outfig = str(OUT_DIR / "kernel_core_test.png")
    fig.savefig(outfig, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {outfig}")
    pl.show()


if __name__ == "__main__":
    main()
