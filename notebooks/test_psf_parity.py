#!/usr/bin/env python3
"""Test script to verify oversampled PSFs are odd-sized and centered."""

import numpy as np
import matplotlib.pyplot as pl
from jwst_kernels.make_psf import save_nircam_PSF, save_miri_PSF, read_PSF
import tempfile
import os


def check_psf_properties(data, label):
    """Print PSF array shape and centering diagnostics."""
    ny, nx = data.shape
    is_odd_x = nx % 2 == 1
    is_odd_y = ny % 2 == 1

    peak_y, peak_x = np.unravel_index(np.argmax(data), data.shape)
    center_y, center_x = (ny - 1) / 2, (nx - 1) / 2

    print(f"\n--- {label} ---")
    print(f"  Shape: {data.shape}")
    print(f"  Odd-sized:  x={is_odd_x}, y={is_odd_y}")
    print(f"  Array center: ({center_y:.1f}, {center_x:.1f})")
    print(f"  Peak pixel:   ({peak_y}, {peak_x})")
    print(f"  Peak offset:  dy={peak_y - center_y:.1f}, dx={peak_x - center_x:.1f}")

    return is_odd_x and is_odd_y, (peak_y == center_y) and (peak_x == center_x)


def main():
    tmpdir = tempfile.mkdtemp(prefix="jwst_psf_test_")
    print(f"Saving PSFs to {tmpdir}\n")

    print("Generating NIRCam F335M PSF...")
    save_nircam_PSF(["F335M"], output_dir=tmpdir)

    print("\nGenerating MIRI F770W PSF...")
    save_miri_PSF(["F770W"], output_dir=tmpdir)

    nircam_filter = {"camera": "NIRCam", "filter": "F335M"}
    miri_filter = {"camera": "MIRI", "filter": "F770W"}

    nircam_data, nircam_pixscale = read_PSF(nircam_filter, psf_dir=tmpdir)
    miri_data, miri_pixscale = read_PSF(miri_filter, psf_dir=tmpdir)

    odd_nircam, centered_nircam = check_psf_properties(nircam_data, "NIRCam F335M (oversampled)")
    odd_miri, centered_miri = check_psf_properties(miri_data, "MIRI F770W (oversampled)")

    print("\n========== RESULTS ==========")
    all_pass = True
    for name, odd, cent in [
        ("NIRCam F335M", odd_nircam, centered_nircam),
        ("MIRI F770W", odd_miri, centered_miri),
    ]:
        status_odd = "PASS" if odd else "FAIL"
        status_cent = "PASS" if cent else "FAIL"
        print(f"  {name}:  odd={status_odd}  centered={status_cent}")
        if not (odd and cent):
            all_pass = False
    print("=============================\n")

    fig, axes = pl.subplots(2, 3, figsize=(14, 9))

    for row, (data, pixscale, label) in enumerate([
        (nircam_data, nircam_pixscale, "NIRCam F335M"),
        (miri_data, miri_pixscale, "MIRI F770W"),
    ]):
        ny, nx = data.shape
        cy, cx = ny // 2, nx // 2

        im = axes[row, 0].imshow(
            np.log10(data / data.max() + 1e-8),
            origin="lower",
            cmap="inferno",
            vmin=-5,
            vmax=0,
        )
        axes[row, 0].set_title(f"{label} — full ({nx}x{ny})")
        axes[row, 0].axhline(cy, color="cyan", ls="--", lw=0.5, alpha=0.7)
        axes[row, 0].axvline(cx, color="cyan", ls="--", lw=0.5, alpha=0.7)
        pl.colorbar(im, ax=axes[row, 0], label="log10(I/Imax)")

        hw = 15
        cutout = data[cy - hw:cy + hw + 1, cx - hw:cx + hw + 1]
        im2 = axes[row, 1].imshow(
            np.log10(cutout / cutout.max() + 1e-8),
            origin="lower",
            cmap="inferno",
            vmin=-3,
            vmax=0,
        )
        axes[row, 1].set_title(f"{label} — center ({2*hw+1}x{2*hw+1})")
        axes[row, 1].axhline(hw, color="cyan", ls="--", lw=0.5, alpha=0.7)
        axes[row, 1].axvline(hw, color="cyan", ls="--", lw=0.5, alpha=0.7)
        pl.colorbar(im2, ax=axes[row, 1], label="log10(I/Imax)")

        row_profile = data[cy, :]
        col_profile = data[:, cx]
        axes[row, 2].plot(np.arange(nx) - cx, row_profile / row_profile.max(), label="row")
        axes[row, 2].plot(np.arange(ny) - cy, col_profile / col_profile.max(), label="col")
        axes[row, 2].set_xlim(-30, 30)
        axes[row, 2].set_yscale("log")
        axes[row, 2].set_ylim(1e-5, 1.5)
        axes[row, 2].axvline(0, color="k", ls=":", lw=0.5)
        axes[row, 2].set_xlabel("Pixel offset from center")
        axes[row, 2].set_ylabel("Normalized intensity")
        axes[row, 2].set_title(f"{label} — cross sections")
        axes[row, 2].legend(fontsize=8)

    fig.suptitle("PSF Parity & Centering Verification", fontsize=14, y=1.01)
    fig.subplots_adjust(wspace=0.4, hspace=0.35)
    tmpdir = '.'
    outpath = os.path.join(tmpdir, "psf_parity_test.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {outpath}")
    pl.show()


if __name__ == "__main__":
    main()
