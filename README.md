# P. falciparum IBD Simulation and Evaluation Framework

A comprehensive forward-time population genetics simulator that explicitly tracks Identity-by-Descent (IBD) segments for evaluating IBD detection methods in *Plasmodium falciparum*.

## 🚀 Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run tests
python3 test_basic_functionality.py

# Run a simple example
python3 -c "from simulation_runner import run_simple_example; run_simple_example()"
```

See [QUICKSTART.md](QUICKSTART.md) for more details.

## 📋 Features

### Core Simulation
- ✅ **Proper chromosome structure** - Realistic 14-chromosome P. falciparum genome
- ✅ **Explicit IBD tracking** - Each segment labeled with founder haplotype
- ✅ **Realistic recombination** - Poisson-distributed crossovers with empirical rates
- ✅ **Multi-population support** - Multiple populations with migration
- ✅ **Demographic events** - Bottlenecks and population size changes
- ✅ **Ground truth IBD** - Perfect knowledge for method evaluation

### IBD Detection Methods
- 🔍 **Threshold-based** - Simple sliding window IBS approach
- 🔍 **HMM-based** - Hidden Markov Model (hmmIBD-style)
- 🔍 **Window-based** - Overlapping windows (isoRelate-style)

### Evaluation Framework
- 📊 Sensitivity, precision, F1 scores
- 📊 Compare across population structures
- 📊 Compare across SNP densities
- 📊 Export results for visualization

## 📁 Project Structure

```
├── ibd_simulation.py        # Core simulation engine
├── ibd_calling.py            # IBD detection methods
├── simulation_runner.py      # High-level interface
├── example_simulations.py    # Example scenarios
├── test_basic_functionality.py  # Unit tests
├── requirements.txt          # Python dependencies
├── README.md                 # This file
├── QUICKSTART.md            # Quick start guide
└── IBD_SIMULATION_README.md # Full documentation
```

## 📖 Documentation

- **[Quick Start Guide](QUICKSTART.md)** - Get started in 5 minutes
- **[Full Documentation](IBD_SIMULATION_README.md)** - Complete API reference and examples

## 🎯 Use Cases

This framework is designed for:

1. **Evaluating IBD detection methods** - Compare algorithm performance with ground truth
2. **Understanding IBD patterns** - How demography affects IBD in malaria parasites
3. **Method development** - Test new IBD calling approaches
4. **Study design** - Determine required SNP density for IBD detection
5. **Population genetics research** - Explore effects of migration, bottlenecks on IBD

## 💡 Example Usage

### Basic Simulation

```python
from simulation_runner import SimulationRunner
import numpy as np

# Create simulator
runner = SimulationRunner(output_dir="my_results")

# Single population of 100 individuals
runner.setup_single_population(pop_size=100)

# Run for 20 generations
individuals = runner.run_simulation(num_generations=20)

# Analyze IBD
sampled = list(np.random.choice(individuals, size=30, replace=False))
true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)
genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)
called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
evaluation = runner.evaluate_methods()
runner.save_results()
```

### Run All Example Scenarios

```bash
python3 example_simulations.py
```

This runs 5 scenarios:
1. Basic single population
2. Population with bottleneck  
3. Two populations with migration
4. Isolated populations
5. SNP density comparison

## 🧪 Testing

```bash
# Run basic functionality tests
python3 test_basic_functionality.py
```

All tests should pass, confirming:
- Segment and chromosome creation
- Population initialization
- Recombination and inheritance
- IBD detection
- Method evaluation

## 📊 Output Files

Each simulation generates:
- `*_individuals.pkl` - Full individual objects
- `*_true_ibd.json` - Ground truth IBD segments
- `*_genotypes_chr*.npy` - Genotype matrices per chromosome
- `*_evaluation.csv` - Detailed evaluation metrics
- `*_summary.csv` - Method summary statistics
- `ibd_proportions.csv` - IBD proportions for plotting

## 🔬 Scientific Background

### P. falciparum Genome
- 14 chromosomes (640 kb - 3.3 Mb)
- Total genome: ~23 Mb
- Recombination rate: ~7.4×10⁻⁷ per bp per generation

### IBD in Malaria
Identity-by-Descent (IBD) segments are genomic regions inherited from a recent common ancestor without recombination. In malaria:
- High IBD indicates recent common ancestry or inbreeding
- Population bottlenecks (epidemics) increase IBD
- IBD patterns reveal transmission dynamics
- Can estimate relatedness and detect outbreaks

## 🛠️ Dependencies

- Python 3.8+
- numpy >= 1.21.0
- pandas >= 1.3.0
- scipy >= 1.7.0
- matplotlib >= 3.4.0 (optional, for plotting)
- seaborn >= 0.11.0 (optional, for plotting)

## 📝 Citation

If you use this simulator, please cite:

```
IBD Segment Tracking Simulator for P. falciparum
Generated for IBD method evaluation
Date: 2025-10-01
```

## 🔗 Related Work

Key papers on IBD detection in malaria:
- Schaffner et al. (2018) - hmmIBD for malaria
- Taylor et al. (2019) - IBD and transmission in P. falciparum
- Machado et al. (2019) - Population structure and IBD
- Hendrickson et al. (2021) - IBD-based relatedness estimation

## 📄 License

This code is provided for research and educational purposes.

## 🤝 Contributing

This is a research tool. Feel free to extend or modify for your needs.

## ⚠️ Note

This is a forward-time simulator designed for moderate-scale simulations (50-200 individuals, 10-50 generations). For very large-scale simulations, consider coalescent-based approaches like msprime.

---

**Created**: October 1, 2025  
**Purpose**: Evaluating IBD detection methods in P. falciparum  
**Status**: Ready for use ✅