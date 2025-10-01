# IBD Segment Tracking Simulator for P. falciparum

A comprehensive forward-time population genetics simulator that explicitly tracks Identity-by-Descent (IBD) segments through recombination and inheritance, designed for evaluating IBD detection methods in *Plasmodium falciparum*.

## Overview

This simulation framework implements:

### Core Features
- **Proper chromosome structure**: Realistic 14-chromosome P. falciparum genome with accurate lengths
- **Explicit IBD tracking**: Each chromosomal segment is labeled with its founder haplotype
- **Recombination**: Poisson-distributed crossovers with realistic recombination rates (~7.4×10⁻⁷ per bp)
- **Multi-population structure**: Support for multiple populations with migration
- **Demographic events**: Bottlenecks and population size changes
- **Ground truth IBD**: Complete knowledge of true IBD segments for method evaluation

### IBD Calling Methods
Three IBD detection methods are implemented for comparison:

1. **Threshold-based caller**: Simple sliding window approach with IBS thresholds
2. **HMM-based caller**: Hidden Markov Model similar to hmmIBD
3. **Window-based caller**: Overlapping windows similar to isoRelate approach

### Evaluation Framework
- Compute sensitivity, precision, and F1 scores against ground truth
- Compare methods across different:
  - Population sizes and structures
  - Demographic histories (bottlenecks)
  - SNP densities
  - Migration rates

## Installation

```bash
# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### Basic Usage

```python
from simulation_runner import SimulationRunner

# Create simulator
runner = SimulationRunner(output_dir="my_results")

# Set up single population
runner.setup_single_population(pop_size=100, recombination_rate=7.4e-7)

# Run simulation for 20 generations
individuals = runner.run_simulation(num_generations=20, allow_selfing=True)

# Sample individuals for analysis
import numpy as np
sampled = list(np.random.choice(individuals, size=30, replace=False))

# Compute ground truth IBD (minimum segment length: 100kb)
true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)

# Generate SNP genotype data
genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)

# Call IBD using all methods
called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)

# Evaluate methods against ground truth
evaluation_df = runner.evaluate_methods()

# Save all results
runner.save_results(prefix="my_simulation")
```

### Run Example Simulations

```bash
# Run all example scenarios
python example_simulations.py

# Or run individual examples in Python
python -c "from example_simulations import example_1_basic_simulation; example_1_basic_simulation()"
```

## Architecture

### File Structure

```
ibd_simulation.py       - Core simulation engine with IBD tracking
ibd_calling.py          - IBD detection methods implementation
simulation_runner.py    - High-level interface for running simulations
example_simulations.py  - Pre-configured example scenarios
requirements.txt        - Python dependencies
```

### Key Classes

#### `ibd_simulation.py`

- **`Segment`**: Represents a chromosomal segment with founder ID
- **`Chromosome`**: A haploid chromosome as a list of IBD segments
- **`Individual`**: Diploid individual with multiple chromosome pairs
- **`ChromosomeStructure`**: P. falciparum chromosome specifications
- **`Population`**: Manages individuals and demographic events
- **`Simulator`**: Multi-population simulation with migration

#### `ibd_calling.py`

- **`ThresholdIBDCaller`**: Simple IBS-based method
- **`HMMIBDCaller`**: Hidden Markov Model approach
- **`WindowIBDCaller`**: Overlapping window method
- **`evaluate_ibd_calling()`**: Performance evaluation functions

#### `simulation_runner.py`

- **`SimulationRunner`**: High-level interface for complete workflows
- Handles setup, execution, calling, evaluation, and output

## Example Scenarios

### 1. Basic Single Population

Simple panmictic population to establish baseline IBD patterns.

```python
from example_simulations import example_1_basic_simulation
example_1_basic_simulation()
```

### 2. Population with Bottleneck

Simulates epidemic-like bottleneck events that increase IBD sharing.

```python
from example_simulations import example_2_with_bottleneck
example_2_with_bottleneck()
```

### 3. Two Populations with Migration

Tests method performance with population structure and gene flow.

```python
from example_simulations import example_3_two_populations_with_migration
example_3_two_populations_with_migration()
```

### 4. Isolated Populations

Creates distinct populations to test detection of population structure.

```python
from example_simulations import example_4_isolated_populations
example_4_isolated_populations()
```

### 5. SNP Density Comparison

Evaluates how SNP density affects IBD calling accuracy.

```python
from example_simulations import example_5_varying_snp_density
example_5_varying_snp_density()
```

## Advanced Usage

### Custom Population Structure

```python
from ibd_simulation import Simulator, ChromosomeStructure
import numpy as np

# Create custom chromosome structure
chrom_structure = ChromosomeStructure(recombination_rate=1e-6)

# Initialize simulator
sim = Simulator(chrom_structure)

# Add multiple populations
sim.add_population(0, size=100)  # Population 0
sim.add_population(1, size=50)   # Population 1
sim.add_population(2, size=75)   # Population 2

# Set migration matrix
migration = np.array([
    [0.90, 0.05, 0.05],  # Pop 0 receives 5% from each other pop
    [0.05, 0.90, 0.05],  # Pop 1
    [0.05, 0.05, 0.90]   # Pop 2
])
sim.set_migration_matrix(migration)

# Define demographic events
events = [
    (10, 'bottleneck', {'population': 0, 'new_size': 20}),
    (15, 'size_change', {'population': 0, 'new_size': 100}),
    (20, 'bottleneck', {'population': 1, 'new_size': 10}),
]

# Run simulation
sim.run_simulation(num_generations=30, demographic_events=events)
```

### Custom IBD Calling Parameters

```python
from ibd_calling import HMMIBDCaller, ThresholdIBDCaller

# Customize HMM parameters
hmm_caller = HMMIBDCaller(
    switch_rate=5e-7,        # IBD state switching rate
    error_rate=0.02,         # Genotyping error rate
    min_segment_length=50000,  # Minimum 50kb segments
    min_posterior=0.7        # Higher confidence threshold
)

# Customize threshold parameters
threshold_caller = ThresholdIBDCaller(
    window_size=100000,      # 100kb windows
    step_size=20000,         # 20kb steps
    ibs_threshold=0.98,      # 98% IBS required
    min_segment_length=200000  # Minimum 200kb segments
)

# Use custom callers
segments = hmm_caller.call_ibd(genotype_data, ind1_idx=0, ind2_idx=1, 
                               chrom_lengths=chrom_structure.lengths)
```

### Access Raw Simulation Data

```python
# Get all individuals from a population
individuals = sim.populations[0].individuals

# Examine chromosome structure
for i, ind in enumerate(individuals[:3]):
    print(f"\nIndividual {ind.id} (Generation {ind.generation}):")
    
    # Look at first chromosome
    maternal, paternal = ind.chromosomes[0]
    
    print(f"  Maternal chromosome: {len(maternal.segments)} segments")
    for seg in maternal.segments[:3]:  # First 3 segments
        print(f"    {seg}")
    
    print(f"  Paternal chromosome: {len(paternal.segments)} segments")
    for seg in paternal.segments[:3]:
        print(f"    {seg}")
```

### Direct IBD Computation

```python
# Compute IBD between specific individuals
ind1 = individuals[0]
ind2 = individuals[5]

# Check maternal haplotypes on chromosome 0
maternal1 = ind1.chromosomes[0][0]
maternal2 = ind2.chromosomes[0][0]

# Find shared IBD regions (minimum 100kb)
ibd_regions = maternal1.get_ibd_with(maternal2, min_length=100000)

print(f"IBD regions between ind {ind1.id} and ind {ind2.id}:")
for start, end in ibd_regions:
    print(f"  {start:.0f} - {end:.0f} bp (length: {(end-start)/1000:.1f} kb)")
```

## Output Files

The simulation generates several output files:

### Per Simulation
- `*_individuals.pkl` - Pickled individual objects (full data)
- `*_true_ibd.json` - Ground truth IBD segments
- `*_genotypes_chr*.npy` - Genotype matrices per chromosome
- `*_evaluation.csv` - Detailed evaluation metrics per pair
- `*_summary.csv` - Summary statistics by method
- `ibd_proportions.csv` - Genome-wide IBD proportions for plotting

### Evaluation Metrics
- **Sensitivity**: Proportion of true IBD correctly identified
- **Precision**: Proportion of called IBD that is truly IBD
- **F1 Score**: Harmonic mean of sensitivity and precision

## P. falciparum Specific Parameters

### Chromosome Lengths (3D7 reference)
```python
CHROMOSOME_LENGTHS = [
    640851,   # Chr 1
    947102,   # Chr 2
    1067971,  # Chr 3
    1200490,  # Chr 4
    1343557,  # Chr 5
    1418242,  # Chr 6
    1445207,  # Chr 7
    1472805,  # Chr 8
    1541735,  # Chr 9
    1687656,  # Chr 10
    2038340,  # Chr 11
    2271494,  # Chr 12
    2925236,  # Chr 13
    3291936,  # Chr 14
]
```

### Recombination Rate
- Default: 7.4×10⁻⁷ per bp per generation
- Based on empirical estimates for P. falciparum
- Can be customized per simulation

### Generation Time
- P. falciparum completes sexual reproduction in mosquito
- One generation ≈ mosquito feeding cycle
- Typical simulations: 10-50 generations

## Design Decisions

### Why Explicit IBD Tracking?

Traditional coalescent simulators work backward in time and may not capture all aspects of IBD accurately in complex demographic scenarios. This forward-time approach:

1. **Directly tracks IBD**: Each segment knows its founder, providing perfect ground truth
2. **Handles complex demography**: Easy to implement bottlenecks, migrations, selection
3. **Realistic recombination**: Breaks up IBD segments naturally over generations
4. **Validation ready**: Ground truth enables rigorous method evaluation

### Computational Considerations

- Forward-time simulation is more computationally intensive than coalescent
- Practical for moderate population sizes (50-200) and generations (10-50)
- Segment representation is memory efficient
- Can be parallelized for multiple replicates

## Future Extensions

Potential additions for enhanced realism:

1. **Within-host dynamics**: Multiple infections, recombination between clones
2. **Selection**: Drug resistance alleles affecting fitness
3. **Non-random mating**: Geographic structure, preferential mating
4. **Variable recombination**: Hotspots and coldspots
5. **Genotyping errors**: Systematic biases, missing data
6. **More IBD callers**: Implement GERMLINE, isoRelate, RefinedIBD equivalents

## Citation

If you use this simulator, please cite:

```
IBD Segment Tracking Simulator for P. falciparum
Author: [Generated for IBD method evaluation]
Date: 2025-10-01
GitHub: [Repository URL]
```

## References

Key papers on IBD detection in malaria:

1. Schaffner et al. (2018) - hmmIBD for malaria
2. Taylor et al. (2019) - IBD and transmission in P. falciparum  
3. Machado et al. (2019) - Population structure and IBD
4. Hendrickson et al. (2021) - IBD-based relatedness estimation

## License

This code is provided for research and educational purposes.

## Contact

For questions or issues, please open a GitHub issue or contact the maintainer.

---

**Last Updated**: October 1, 2025
