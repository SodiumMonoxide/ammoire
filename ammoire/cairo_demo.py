"""Export the Cairo unstrained/strained bands and optional splitting maps."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .models import make_model


def spin_bands(model, k):
    h = model.hamiltonian(k)
    return np.stack([np.linalg.eigvalsh(h[..., s::2, s::2]) for s in (0, 1)], axis=-2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strain', type=float, default=.1, help='diagonal hopping anisotropy delta')
    parser.add_argument('--points', type=int, default=60, help='points per path segment')
    parser.add_argument('--out', type=Path, default=Path('results/cairo'))
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.points < 2 or not np.isfinite(args.strain):
        parser.error('points must be >=2 and strain finite')
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            parser.error('--plot needs matplotlib (install ammoire[plot])')
    vertices = np.pi * np.array([[0, 0], [1, 0], [1, 1], [0, 0]])
    points = np.vstack([np.linspace(a, b, args.points, endpoint=False)
                        for a, b in zip(vertices[:-1], vertices[1:])] + [vertices[-1:]])
    args.out.mkdir(parents=True, exist_ok=True)
    cases = [('unstrained', make_model('cairo')), ('strained', make_model('cairo', delta=args.strain))]
    summary = {'source': cases[0][1].source, 'unit': '|t| (paper units)',
               'basis': 'A,B,alpha,beta,gamma,eta; physical spin interleaved',
               'gauge': 'periodic cell gauge, d=0 source transformed by U',
               'strained_delta': args.strain, 'band_comparison': 'within-spin energy rank, not wavefunction tracking'}
    with (args.out / 'bands.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['case', 'kx_a', 'ky_a', 'ordered_band_1based', 'up', 'down', 'up_minus_down'])
        for label, model in cases:
            bands = spin_bands(model, points)
            for k, spectrum in zip(points, bands):
                for band, (up, down) in enumerate(spectrum.T, 1):
                    writer.writerow([label, *k, band, up, down, up - down])
            r = np.array([.04, .02])
            near = spin_bands(model, np.stack([2 * r, r], axis=-1))
            delta = near[:, 0, 4] - near[:, 1, 4]
            power = float(np.log2(abs(delta[0] / delta[1]))) if np.all(abs(delta) > 1e-14) else None
            summary[label] = {'fifth_branch_gamma_power': power,
                              'radius_values_ka': r.tolist(), 'splitting_on_2r_r': delta.tolist()}
    if args.plot:
        axis = np.linspace(-np.pi, np.pi, 101)
        x, y = np.meshgrid(axis, axis)
        k = np.stack([x, y], axis=-1)
        maps = []
        for _, model in cases:
            e = spin_bands(model, k)
            maps.append(e[..., 0, 4] - e[..., 1, 4])
        limit = max(np.max(abs(field)) for field in maps)
        fig, panels = plt.subplots(1, 2, figsize=(9, 4), layout='constrained')
        for panel, field, title in zip(panels, maps, ['Unstrained', f'Diagonal strain delta={args.strain:g}']):
            im = panel.imshow(field, origin='lower', extent=[-1, 1, -1, 1],
                              cmap='RdBu_r', vmin=-limit, vmax=limit)
            panel.set(xlabel='kx / pi', ylabel='ky / pi', title=title)
        fig.colorbar(im, ax=panels, label='Fifth ordered branch: up - down (|t|)')
        fig.savefig(args.out / 'splitting.png', dpi=160)
        plt.close(fig)
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
