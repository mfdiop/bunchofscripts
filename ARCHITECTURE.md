# Architecture Overview

## System Design

The IBD simulation framework consists of four main modules:

```
┌─────────────────────────────────────────────────────────────┐
│                     User Interface Layer                     │
│                                                               │
│  ┌────────────────────┐        ┌─────────────────────────┐  │
│  │ example_           │        │ Custom Scripts          │  │
│  │ simulations.py     │        │ (User Code)             │  │
│  └────────┬───────────┘        └───────────┬─────────────┘  │
│           │                                 │                │
│           └────────────────┬────────────────┘                │
└────────────────────────────┼─────────────────────────────────┘
                             │
┌────────────────────────────┼─────────────────────────────────┐
│                High-Level Interface Layer                     │
│                                                               │
│           ┌────────────────┴──────────────────┐              │
│           │   simulation_runner.py             │              │
│           │                                    │              │
│           │  - SimulationRunner class          │              │
│           │  - Workflow management             │              │
│           │  - Output formatting               │              │
│           └──────────┬──────────────┬──────────┘              │
└──────────────────────┼──────────────┼───────────────────────┘
                       │              │
        ┌──────────────┘              └──────────────┐
        │                                            │
┌───────┼────────────────────┐        ┌──────────────┼─────────┐
│  Core Simulation Layer     │        │  IBD Calling Layer     │
│                            │        │                        │
│  ibd_simulation.py         │        │  ibd_calling.py        │
│                            │        │                        │
│  ┌──────────────────┐      │        │  ┌──────────────────┐ │
│  │ Segment          │      │        │  │ ThresholdCaller  │ │
│  │ Chromosome       │      │        │  │ HMMIBDCaller     │ │
│  │ Individual       │      │        │  │ WindowCaller     │ │
│  └──────────────────┘      │        │  └──────────────────┘ │
│                            │        │                        │
│  ┌──────────────────┐      │        │  ┌──────────────────┐ │
│  │ Population       │      │        │  │ evaluate_ibd     │ │
│  │ Simulator        │      │        │  │ _calling()       │ │
│  └──────────────────┘      │        │  └──────────────────┘ │
│                            │        │                        │
│  ┌──────────────────┐      │        └────────────────────────┘
│  │ ChromosomeStruct │      │
│  └──────────────────┘      │
└─────────────────────────────┘
```

## Data Flow

### 1. Simulation Phase

```
Initialize Founders
       │
       ├─> Each founder gets unique haplotype labels
       │
       v
  Generation Loop
       │
       ├─> For each individual in next generation:
       │   │
       │   ├─> Select two parents
       │   │
       │   ├─> Meiosis 1: Parent 1
       │   │   ├─> Generate recombination breakpoints
       │   │   └─> Create gamete with IBD segments
       │   │
       │   ├─> Meiosis 2: Parent 2
       │   │   ├─> Generate recombination breakpoints
       │   │   └─> Create gamete with IBD segments
       │   │
       │   └─> Combine gametes → New individual
       │
       ├─> Apply demographic events (bottlenecks, etc.)
       │
       └─> Continue for N generations

Result: Individuals with explicit IBD structure
```

### 2. Genotyping Phase

```
True IBD Structure
       │
       ├─> Place SNPs along chromosomes
       │
       ├─> For each individual:
       │   └─> For each SNP position:
       │       └─> Determine genotype from founder allele
       │
       v
Genotype Matrix
(SNPs × Individuals)
```

### 3. IBD Calling Phase

```
Genotype Data
       │
       ├──────────┬──────────┬──────────┐
       │          │          │          │
       v          v          v          v
  Threshold    HMM     Window    [Custom]
   Method     Method   Method     Method
       │          │          │          │
       └──────────┴──────────┴──────────┘
                     │
                     v
           Called IBD Segments
```

### 4. Evaluation Phase

```
Called IBD ──┐
             ├──> Compare Segments ──> Metrics
True IBD ────┘         │             (Sensitivity,
                       │              Precision,
                       v              F1 score)
              Per-base comparison
                       │
                       ├─> True Positives
                       ├─> False Positives
                       └─> False Negatives
```

## Key Data Structures

### Segment
```python
class Segment:
    start: float           # Position in bp
    end: float            # Position in bp
    founder_id: int       # Original founder haplotype
    generation_origin: int # When this segment was created
```

### Chromosome
```python
class Chromosome:
    segments: List[Segment]  # Non-overlapping, contiguous
    length: float           # Total chromosome length
```

**Invariants:**
- Segments are sorted by position
- No gaps between segments
- No overlapping segments
- Segments cover entire chromosome [0, length]

### Individual
```python
class Individual:
    id: int
    chromosomes: List[Tuple[Chromosome, Chromosome]]  # [maternal, paternal]
    generation: int
    population: int
```

**Structure:**
- 14 chromosome pairs (P. falciparum)
- Each pair: (maternal, paternal)
- Each chromosome: list of IBD segments

## Recombination Algorithm

```
Input: Two parental chromosomes, breakpoint positions

1. Initialize empty segment list
2. Set current_parent = 0 (first parent)
3. For each region between breakpoints:
   a. Get segments from current parent in this region
   b. For each segment:
      - Find overlap with current region
      - Add overlapping portion to new segment list
   c. Switch current_parent (0 ↔ 1)
4. Merge adjacent segments with same founder_id
5. Return new chromosome

Output: Recombinant chromosome with IBD tracking
```

**Example:**
```
Parent 1: [═══A═══|═══B═══|═══C═══]
Parent 2: [───X───|───Y───|───Z───]
                  ↑       ↑
           Breakpoints

Result:   [═══A═══|───Y───|═══C═══]
          (A from P1, Y from P2, C from P1)
```

## IBD Detection: HMM Approach

```
States: {IBD, not-IBD}

Emission Probabilities:
- P(match | IBD) = 1 - error_rate
- P(match | not-IBD) = p² + (1-p)²  (by chance)
  where p = allele frequency

Transition Probabilities:
- P(switch) = 1 - exp(-λ × distance)
  where λ = switch_rate

Algorithm:
1. Compute emission probabilities at each SNP
2. Run forward algorithm (α)
3. Run backward algorithm (β)  
4. Compute posterior: γ = α × β (normalized)
5. Call IBD where γ(IBD) > threshold
```

## Performance Characteristics

### Time Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Initialize founders | O(N × C) | N=pop size, C=14 chromosomes |
| One mating | O(C × R × S) | R=recombination events, S=segments |
| One generation | O(N × C × R × S) | N matings |
| IBD computation | O(N² × C × S) | Pairwise, all chromosomes |
| Threshold calling | O(M × W) | M=SNPs, W=windows |
| HMM calling | O(M × K) | K=states (2) |

### Space Complexity

| Structure | Complexity | Notes |
|-----------|-----------|-------|
| One individual | O(C × S) | C=14, S=avg segments |
| Population | O(N × C × S) | N individuals |
| Genotype data | O(M × N) | M SNPs, N individuals |
| IBD results | O(N² × I) | I=IBD segments per pair |

### Practical Limits

**Recommended:**
- Population size: 50-200
- Generations: 10-50
- SNPs: 1000-5000 per chromosome
- Sample size: 20-100 for analysis

**Maximum (single machine):**
- Population size: ~500
- Generations: ~100
- SNPs: ~10,000 per chromosome

## Extension Points

### Adding New IBD Calling Methods

1. Create new class in `ibd_calling.py`
2. Implement `call_ibd()` method
3. Return list of `CalledIBDSegment` objects
4. Add to `simulation_runner.py` caller dictionary

```python
class MyCustomCaller:
    def call_ibd(self, genotype_data, ind1_idx, ind2_idx, 
                 chrom_lengths) -> List[CalledIBDSegment]:
        # Your algorithm here
        return segments
```

### Adding Demographic Events

Add to `Population` class in `ibd_simulation.py`:

```python
def apply_custom_event(self, params):
    # Modify self.individuals
    # Update self.size if needed
    pass
```

Register in `Simulator.run_simulation()`:
```python
if event_type == 'custom_event':
    pop.apply_custom_event(params)
```

### Adding Selection

Modify `Population.advance_generation()`:

```python
# Assign fitness based on genotype
fitnesses = self.compute_fitness(self.individuals)

# Sample parents with probability ∝ fitness
parents = np.random.choice(
    self.individuals, 
    p=fitnesses/fitnesses.sum()
)
```

### Adding Complex Mating

Modify `Simulator.advance_generation_with_migration()`:

```python
# Custom parent selection
parent1, parent2 = self.select_parents_custom(
    population, 
    geographic_structure,
    mating_preferences
)
```

## Testing Strategy

### Unit Tests (`test_basic_functionality.py`)
- Segment operations
- Chromosome creation and IBD detection
- Population initialization
- Recombination mechanics
- IBD calling methods

### Integration Tests (Example simulations)
- End-to-end workflows
- Method comparisons
- Output validation

### Validation Tests
- Known IBD patterns (e.g., parent-offspring should be ~50% IBD)
- Migration equilibrium (F_ST calculations)
- Bottleneck effects (increased IBD)

## Future Enhancements

1. **Performance**
   - Cython implementation of hot loops
   - Parallel population simulations
   - Vectorized recombination

2. **Features**
   - Within-host mixed infections
   - Drug resistance mutations
   - Geographic structure
   - Seasonal fluctuations

3. **Output**
   - VCF format export
   - Visualization functions
   - Interactive dashboards

4. **Validation**
   - Compare to real P. falciparum data
   - Benchmark against coalescent simulators
   - Statistical tests for accuracy
