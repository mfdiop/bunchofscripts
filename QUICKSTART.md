# Quick Start Guide

Get up and running with the IBD simulator in 5 minutes.

## Installation

```bash
# Install dependencies
pip install -r requirements.txt
```

## Run Tests

```bash
# Verify installation
python3 test_basic_functionality.py
```

## Run Your First Simulation

### Option 1: Simple Example (Fastest)

```bash
python3 -c "from simulation_runner import run_simple_example; run_simple_example()"
```

This will:
- Create a population of 50 individuals
- Simulate for 10 generations
- Sample 20 individuals
- Compute ground truth IBD
- Call IBD using 3 methods (Threshold, HMM, Window)
- Evaluate method performance
- Save results to `example_results/`

### Option 2: Run Example Scenarios

```bash
# Run all example scenarios
python3 example_simulations.py
```

This runs 5 different scenarios:
1. Basic single population
2. Population with bottleneck
3. Two populations with migration
4. Isolated populations
5. SNP density comparison

Results saved to separate directories (`results_example1/`, etc.)

### Option 3: Custom Simulation

```python
from simulation_runner import SimulationRunner
import numpy as np

# Initialize
runner = SimulationRunner(output_dir="my_simulation")

# Set up population
runner.setup_single_population(pop_size=100)

# Run for 20 generations
individuals = runner.run_simulation(num_generations=20, allow_selfing=True)

# Sample some individuals
sampled = list(np.random.choice(individuals, size=30, replace=False))

# Compute ground truth
true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)

# Generate genotype data
genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)

# Call IBD with all methods
called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)

# Evaluate
evaluation_df = runner.evaluate_methods()
print(evaluation_df.groupby('method')[['sensitivity', 'precision', 'f1_score']].mean())

# Save results
runner.save_results(prefix="my_sim")
```

## Understanding Output

### Key Output Files

```
results_directory/
├── *_individuals.pkl        # Full individual objects
├── *_true_ibd.json          # Ground truth IBD segments
├── *_genotypes_chr*.npy     # Genotype matrices
├── *_evaluation.csv         # Detailed per-pair metrics
├── *_summary.csv            # Method summary statistics
└── ibd_proportions.csv      # IBD proportions for plotting
```

### Reading Results

```python
import pandas as pd

# Read evaluation results
eval_df = pd.read_csv("results_example1/basic_simulation_evaluation.csv")

# Summary by method
summary = eval_df.groupby('method')[['sensitivity', 'precision', 'f1_score']].agg(['mean', 'std'])
print(summary)

# Read IBD proportions
ibd_df = pd.read_csv("results_example1/ibd_proportions.csv")
print(f"Average IBD proportion: {ibd_df['ibd_proportion'].mean():.4f}")
```

## Common Use Cases

### 1. Test Different Population Sizes

```python
runner = SimulationRunner()

for size in [50, 100, 200]:
    runner.setup_single_population(pop_size=size)
    individuals = runner.run_simulation(num_generations=15)
    # ... analyze results
```

### 2. Simulate Bottleneck

```python
runner = SimulationRunner()
runner.setup_single_population(pop_size=200)

# Bottleneck at generation 10
events = [
    (10, 'bottleneck', {'population': 0, 'new_size': 20}),
    (11, 'size_change', {'population': 0, 'new_size': 100})
]

individuals = runner.run_simulation(num_generations=30, demographic_events=events)
```

### 3. Multi-Population with Migration

```python
import numpy as np

runner = SimulationRunner()

# Two populations with 5% migration
migration_matrix = np.array([
    [0.95, 0.05],
    [0.05, 0.95]
])

runner.setup_multi_population(
    pop_sizes=[80, 60],
    migration_matrix=migration_matrix
)

individuals = runner.run_simulation(num_generations=25)
```

### 4. Compare SNP Densities

```python
runner = SimulationRunner()
runner.setup_single_population(pop_size=100)
individuals = runner.run_simulation(num_generations=15)
sampled = individuals[:30]

for snp_count in [100, 500, 1000, 2000]:
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=snp_count)
    called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
    # ... compare results
```

## Next Steps

- Read the [full documentation](IBD_SIMULATION_README.md)
- Explore `example_simulations.py` for detailed examples
- Modify parameters to match your research questions
- Implement custom IBD calling methods in `ibd_calling.py`

## Troubleshooting

### ImportError: No module named 'numpy'

```bash
pip install -r requirements.txt
```

### Simulation runs slowly

- Reduce population size (default: 100 → try 50)
- Reduce number of generations (default: 20 → try 10)
- Sample fewer individuals for analysis

### Low IBD detection

- Increase number of generations to build up IBD
- Add a bottleneck event to increase IBD sharing
- Adjust `min_length` parameter when calling IBD

### Method evaluation gives unexpected results

- Check SNP density (need enough SNPs for detection)
- Verify `min_segment_length` is appropriate
- Examine raw IBD calls vs ground truth in detail

## Support

For issues or questions:
1. Check the [full documentation](IBD_SIMULATION_README.md)
2. Review example code in `example_simulations.py`
3. Run `test_basic_functionality.py` to verify installation
