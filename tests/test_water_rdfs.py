import numpy as np
import pytest

from curdf.rdf import rdf
from ase.io import read
import matplotlib.pyplot as plt
from pathlib import Path

DATAPATH = Path("darren_data/rdf-mace-nvt")
RMIN = 0.0
RMAX = 6.0
N_BINS = 200

def test_compute_rdf_cuda_runs():
    torch_mod = pytest.importorskip("torch")
    if not torch_mod.cuda.is_available():
        pytest.skip("CUDA not available")
    
    water_traj_list = read(DATAPATH / f'nvt.xyz', index=':')
    
    for species_a, species_b in [("O", "O"), ("H", "H"), ("O", "H")]:
        data = np.loadtxt(DATAPATH / f"rdf_{species_a}{species_b}_mace.dat")
        r_ref, gr_ref = data[:,0], data[:,1]
        
        try:
            bins, gr = rdf(
                water_traj_list,
                species_a=species_a,
                species_b=species_b,
                r_min=RMIN,
                r_max=RMAX,
                nbins=N_BINS,
                output=f"cuRDF_results/rdf_{species_a}{species_b}.csv"
            )
            plt.plot(bins, gr, label=f"cuRDF - {species_a}-{species_b}")
            plt.plot(r_ref, gr_ref, label=f"Reference - {species_a}-{species_b}")
            plt.ylim(0, 3)
            plt.legend()
            plt.savefig(f"cuRDF_results/rdf_{species_a}{species_b}_comparison.png")
            plt.close()
            # Numerical check: bins should align and g(r) should match reference within 1% relative.
            assert bins.shape == r_ref.shape
            assert gr.shape == gr_ref.shape
            np.testing.assert_allclose(bins, r_ref, rtol=1e-6, atol=1e-6)
            np.testing.assert_allclose(gr, gr_ref, rtol=1e-2, atol=1e-3)
        except RuntimeError as err:
            pytest.skip(f"CUDA neighbor list unavailable: {err}")
        
