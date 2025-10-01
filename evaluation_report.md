# Detailed Evaluation of Three P. falciparum Simulation Functions

## Overview
This document provides a comprehensive evaluation of three R functions designed to simulate *Plasmodium falciparum* genomic data with the goal of creating ground truth for evaluating genetic relatedness and transmission dynamics methods.

---

## Code 1: `simulate_pf()` - Generation-Based Transmission Model

### Strengths

1. **Explicit Transmission Chain**: Tracks parent-offspring relationships through the `parent` field in metadata, creating a genealogy that can be used as ground truth for transmission inference.

2. **Biologically Realistic COI Distribution**: Uses empirically-derived probabilities for COI (1-6) with realistic weights (55% monoclonal, decreasing for higher COI).

3. **Sexual Recombination**: Implements mating between two parasites (donor and mate) before generating offspring, mimicking the mosquito stage of the parasite lifecycle.

4. **Population Evolution**: Maintains a growing population over time, with new infections arising from existing ones, allowing for realistic population structure to develop.

5. **Multiple Output Formats**: Provides both individual-level data (host allele frequencies) and population-level summaries (genotype matrix, pairwise distances).

### Weaknesses

1. **Unrealistic Recombination Model**:
   - Uses Poisson-distributed crossover count based on entire genome length (1 Mb)
   - For *P. falciparum* (14 chromosomes, ~23 Mb total), recombination occurs during meiosis with ~1-2 crossovers per chromosome
   - Current model doesn't account for chromosome structure or realistic crossover rates

2. **Mutation Rate Application**:
   - Applies mutation uniformly per SNP per generation
   - Doesn't distinguish between synonymous/non-synonymous mutations or account for selection pressure
   - For SNPs, mutation rate of 1e-7 per base per generation is reasonable but ignores microsatellite hypervariability

3. **Overly Simplified Transmission**:
   - No transmission bottleneck (all clones from donor can potentially contribute)
   - Ignores the complex mosquito-stage dynamics where only gametocytes participate
   - No consideration of geographic/household clustering
   - No immigration/migration from outside population

4. **COI Generation**:
   - New infection COI is independent of donor/mate COI
   - In reality, COI relates to superinfection, co-transmission, and local transmission intensity

5. **Consensus Genotype Issues**:
   - Mixed infections (COI > 1) result in many NA values when alleles differ
   - This drastically reduces information content for downstream analyses
   - Doesn't use within-host allele frequencies (WSAF) which is standard for polyclonal samples

6. **Pairwise Hamming Distance**:
   - Simple distance metric that doesn't account for:
     - Population allele frequencies (no correction for IBS vs IBD)
     - Missing data from mixed infections
     - Recombination patterns that would better reflect relatedness

7. **No Migration or Population Structure**:
   - Assumes closed population (island model mentioned in request)
   - No spatial structure or subpopulations

8. **Missing Key Evolutionary Forces**:
   - No selection (all alleles evolve neutrally)
   - No linkage disequilibrium modeling beyond what naturally arises
   - No variation in transmission intensity over time

### Biological Realism Assessment

| Feature | Present | Realistic | Comments |
|---------|---------|-----------|----------|
| COI distribution | ✓ | ✓ | Good empirical distribution |
| Recombination | ✓ | ✗ | Incorrect model for crossovers |
| Mutation | ✓ | ~ | Rate reasonable but no selection |
| Mixed infections | ✓ | ✓ | COI handled correctly |
| Transmission chain | ✓ | ✗ | Too simplified, no bottleneck |
| Migration | ✗ | ✗ | Absent |
| Selection | ✗ | ✗ | Absent |
| LD structure | ~ | ~ | Emerges but not explicitly modeled |
| Transmission intensity | ✗ | ✗ | Fixed over time |

---

## Code 2: `simulate_population()` - Static Population Model

### Strengths

1. **Flexible COI Distribution**: Uses negative binomial distribution with dispersion parameter, allowing control over overdispersion commonly seen in malaria COI data.

2. **Within-Host Allele Frequencies (WSAF)**: Computes allele frequencies across clones within each host, which is the standard representation for polyclonal infections.

3. **VCF-like Output Format**: Generates GT:AD:DP format strings mimicking standard VCF files, making output compatible with existing analysis pipelines.

4. **IBD Ground Truth**: Attempts to compute a ground truth IBD matrix based on shared haplotypes.

5. **Cleaner Structure**: More concise code that's easier to understand and modify.

### Weaknesses

1. **No Temporal Structure**:
   - Single time point snapshot
   - Cannot evaluate methods that use temporal patterns
   - No epidemiological dynamics

2. **Static Founder Pool**:
   - All infections sample from same 50 founder haplotypes
   - No evolution or accumulation of diversity over time
   - Unrealistic for 6-month period mentioned in parameters (though not actually used)

3. **Oversimplified Recombination**:
   - Fixed lambda=2 crossovers regardless of genome size
   - All SNPs treated as if on single chromosome
   - Breakpoints chosen randomly across SNP positions without regard to physical distance
   - No interference or realistic crossover distribution

4. **Flawed IBD Calculation**:
   - Compares haplotypes as concatenated strings (exact matches only)
   - In reality, IBD is measured as proportion of genome shared through recent common ancestry
   - Method won't detect partial sharing or segment-based IBD
   - Division by `max(coi[i], coi[j])` is arbitrary and not biologically motivated

5. **Fixed Coverage**:
   - All SNPs have exactly 20x coverage
   - Real data has variable coverage following negative binomial or Poisson distributions
   - No missing data or low-quality calls

6. **Simplistic Genotype Calling**:
   - Hard thresholds at 0.1 and 0.9 for homozygous calls
   - Real genotype calling uses probabilistic methods accounting for error rates
   - No heterozygosity consideration for diploid organisms (though P. falciparum is haploid)

7. **No Population Structure**:
   - All hosts sample randomly from same pool
   - No relatedness structure, household clustering, or geographic patterns
   - Makes IBD ground truth questionable as measure of recent transmission

8. **Missing Evolutionary Forces**:
   - No mutation
   - No selection
   - No migration
   - Linkage disequilibrium only from finite sampling, not biological processes

### Biological Realism Assessment

| Feature | Present | Realistic | Comments |
|---------|---------|-----------|----------|
| COI distribution | ✓ | ✓ | Negative binomial is appropriate |
| Recombination | ✓ | ✗ | Too simplified, no chromosome structure |
| Mutation | ✗ | ✗ | Completely absent |
| Mixed infections | ✓ | ✓ | WSAF representation good |
| Transmission chain | ✗ | ✗ | No transmission modeled |
| Migration | ✗ | ✗ | Absent |
| Selection | ✗ | ✗ | Absent |
| LD structure | ✗ | ✗ | Not modeled |
| Transmission intensity | ✗ | ✗ | Not applicable |

---

## Code 3: `simulate_population_time()` - Daily Sampling Model

### Strengths

1. **Temporal Resolution**: Daily time steps with stochastic sampling, providing realistic longitudinal data structure.

2. **Flexible Sampling**: Poisson-distributed daily samples allows for realistic variability in sampling effort.

3. **Monthly Aggregation**: Summarizes IBD by time periods, useful for evaluating temporal patterns in relatedness.

4. **Sample Metadata**: Tracks sampling day and derived month for each sample, enabling time-stratified analyses.

5. **VCF-like Output**: Maintains GT:AD:DP format from Code 2.

### Weaknesses

1. **Identical Core Issues to Code 2**:
   - Same flawed recombination model
   - Same problematic IBD calculation
   - Same static founder pool
   - No mutation, selection, or migration

2. **No Transmission Dynamics**:
   - Despite temporal structure, samples are independent draws from same founder pool
   - No epidemiological model linking samples over time
   - Cannot generate realistic transmission chains or outbreak patterns
   - Missing the entire point of temporal data: seeing evolution and transmission links

3. **Unrealistic Evolution**:
   - Population doesn't evolve over 180 days
   - All samples equally related to founders regardless of sampling time
   - No accumulation of diversity
   - No seasonal or temporal trends in transmission intensity

4. **IBD Summary Limitations**:
   - Uses `month_pair = pmax(month_i, month_j)` which aggregates by the later month only
   - Doesn't distinguish within-month vs between-month comparisons
   - Misses opportunity to show decay of relatedness over time

5. **Computational Inefficiency**:
   - Computes all pairwise IBD even for large samples (O(n²) complexity)
   - Stores full genotype and IBD matrices in memory

6. **Same Missing Elements**:
   - No mutation accumulation over time
   - No selection pressure
   - No population structure changes
   - No migration events
   - No variation in transmission intensity

### Biological Realism Assessment

| Feature | Present | Realistic | Comments |
|---------|---------|-----------|----------|
| COI distribution | ✓ | ✓ | Same good NB distribution |
| Recombination | ✓ | ✗ | Same issues as Code 2 |
| Mutation | ✗ | ✗ | Critical absence for temporal data |
| Mixed infections | ✓ | ✓ | WSAF representation good |
| Transmission chain | ✗ | ✗ | Completely missing despite temporal structure |
| Migration | ✗ | ✗ | Absent |
| Selection | ✗ | ✗ | Absent |
| LD structure | ✗ | ✗ | Not modeled |
| Transmission intensity | ✗ | ✗ | Fixed over time |

---

## Comparative Summary

| Aspect | Code 1 | Code 2 | Code 3 |
|--------|--------|--------|--------|
| **Best for** | Transmission chain inference | Population genetic structure | Temporal trend analysis |
| **Transmission model** | Basic genealogy | None | None |
| **Temporal structure** | Monthly bins | Single snapshot | Daily sampling |
| **Output format** | Genotype matrix + metadata | VCF-like | VCF-like |
| **IBD ground truth** | Hamming distance | Shared haplotypes | Shared haplotypes |
| **Population dynamics** | Growing population | Static | Static |
| **Code complexity** | High | Low | Medium |
| **Computational cost** | Medium | Low | Medium-High |
| **Biological realism** | Low-Medium | Low | Low |

---

## Critical Gaps Across All Three Functions

### 1. **Recombination Modeling**
- **Current**: Simplified crossover models without chromosome structure
- **Needed**: 
  - 14 chromosomes with independent segregation
  - Realistic crossover rates per chromosome (~1-2 per meiosis)
  - No crossovers in subtelomeric regions
  - Proper genetic distance vs physical distance mapping

### 2. **Transmission Bottlenecks**
- **Current**: Unrestricted sampling from donor parasites
- **Needed**:
  - Mosquito stage bottleneck (only gametocytes transmit)
  - Limited number of oocysts per mosquito (~1-10)
  - Founding population size in new host
  - Clonal vs mixed transmission events

### 3. **Selection and Drug Resistance**
- **Current**: All loci evolve neutrally
- **Needed**:
  - Known resistance loci (pfcrt, pfmdr1, pfdhfr, pfdhps, pfkelch13)
  - Selection coefficients based on drug pressure
  - Epistatic interactions
  - Fitness costs of resistance

### 4. **Migration/Gene Flow**
- **Current**: Closed island population
- **Needed**:
  - Immigration rate from external reservoir
  - Spatial structure with local transmission clusters
  - Household/village clustering
  - Distance-decay relationships

### 5. **Linkage Disequilibrium**
- **Current**: Emerges passively or absent
- **Needed**:
  - Initial LD structure reflecting population history
  - Decay with recombination
  - LD blocks around selection targets
  - Validation against real *P. falciparum* LD patterns

### 6. **Transmission Intensity Variation**
- **Current**: Fixed parameters throughout simulation
- **Needed**:
  - Seasonal variation in transmission
  - Intervention effects (bed nets, drugs, vaccines)
  - Stochastic variation in local transmission
  - Feedback between prevalence and force of infection

### 7. **IBD Calculation**
- **Current**: String matching (Code 2,3) or Hamming distance (Code 1)
- **Needed**:
  - Segment-based IBD calling
  - Hidden Markov Model approaches
  - Minimum segment length thresholds
  - Proper handling of mixed infections
  - Methods comparable to hmmIBD, isoRelate, or CHIIMP

### 8. **Realistic Genomic Features**
- **Current**: Random SNPs with uniform properties
- **Needed**:
  - Real *P. falciparum* genome structure (23 Mb, 14 chromosomes)
  - Heterogeneous mutation rates (GC content effects)
  - Centromeres and telomeres
  - Regions of low complexity (microsatellites, var genes)
  - Realistic minor allele frequency spectrum

### 9. **Superinfection Dynamics**
- **Current**: Each sample independent
- **Needed**:
  - Within-host competition between clones
  - Sequential infections in same host
  - Immune clearance of genetically similar parasites
  - Duration of infection and chronicity

### 10. **Validation Against Real Data**
- **Current**: No validation
- **Needed**:
  - Compare simulated vs real data for:
    - Allele frequency spectra
    - LD decay curves
    - COI distributions
    - Genetic diversity metrics (π, Tajima's D)
    - IBD sharing distributions
    - Spatial autocorrelation

---

## Recommendations for Improvement

### Priority 1: Essential for Ground Truth IBD

1. **Implement proper IBD tracking**: Record actual genomic segments inherited IBD during simulation, rather than inferring post-hoc from genotypes.

2. **Add realistic recombination**: Model 14 chromosomes with appropriate crossover distributions.

3. **Create transmission chains**: Explicitly model who-infected-whom with timing information.

4. **Add mutation over time**: Essential for temporal data to accumulate divergence.

### Priority 2: Important for Biological Realism

5. **Model transmission bottlenecks**: Limit clone diversity in transmission events.

6. **Add migration**: Import external genetic diversity to prevent unrealistic population structure.

7. **Variable coverage**: Simulate realistic sequencing depth variation and missing data.

8. **Spatial structure**: Add household/village clustering of cases.

### Priority 3: Advanced Features

9. **Selection at key loci**: Model drug resistance evolution.

10. **Seasonal transmission**: Vary intensity over time.

11. **Superinfection**: Allow hosts to accumulate multiple infections.

12. **Validation metrics**: Compare systematically to real *P. falciparum* data.

---

## Suitability for Stated Goals

### Goal: Simulate ground truth relatedness/IBD and evaluate methods

**Code 1**: 
- ✓ Has transmission chains (but too simplified)
- ✗ IBD calculation is just Hamming distance
- ✗ Missing many key biological features
- **Verdict**: Could work for basic method comparison but needs substantial improvement

**Code 2**: 
- ✗ No transmission structure
- ✗ Flawed IBD calculation
- ✗ Static population unrealistic
- **Verdict**: Not suitable without major revision

**Code 3**: 
- ✓ Has temporal structure (but no actual evolution)
- ✗ Same flawed IBD calculation
- ✗ No transmission chains despite temporal data
- **Verdict**: Missed opportunity; temporal structure wasted without transmission model

### Overall Assessment

**None of the three functions are currently suitable** for rigorous evaluation of IBD/relatedness methods because:

1. **Ground truth IBD is not properly defined or tracked** - should be recorded during simulation, not inferred post-hoc
2. **Transmission dynamics are absent or oversimplified** - cannot generate realistic patterns of relatedness
3. **Critical evolutionary forces missing** - mutation, selection, migration all absent or inadequate
4. **Recombination models unrealistic** - will produce wrong IBD segment length distributions

---

## Recommended Next Steps

### Option A: Extend Code 1 (Moderate Effort)
- Add proper chromosome structure to recombination
- Track true IBD segments during inheritance
- Add migration and bottlenecks
- Implement proper IBD calling methods to compare against ground truth

### Option B: Start Fresh with Forward-Time Simulator (High Effort)
- Use existing framework like SLiM or MSprime for core evolution
- Build *P. falciparum*-specific wrapper for:
  - Multi-chromosome genome
  - Complex lifecycle (mosquito/human)
  - Mixed infections and COI
  - Spatial structure
- Record ground truth IBD segments
- Generate realistic sequencing data

### Option C: Use/Extend Existing Malaria Simulators (Low-Medium Effort)
- Leverage existing tools:
  - **OpenMalaria**: Individual-based model with genetics
  - **MOIRE**: Relatedness inference with simulation framework
  - **Coil**: COI estimation with simulation tools
  - **SimRelatedness**: R package for simulating related samples
- Extend to output format needed for evaluation
- Validate against real data

---

## Conclusion

While all three functions demonstrate understanding of key concepts (COI, recombination, mixed infections), **none are adequate** for rigorous method evaluation in their current form. The fundamental issue is that **ground truth IBD must be explicitly tracked during simulation**, not inferred afterward using the same methods being evaluated.

For your stated goal of evaluating methods for determining genetic relatedness in *P. falciparum*, I recommend:

1. **Immediate**: Use Code 1 as a starting point but with major revisions (track IBD segments explicitly)
2. **Short-term**: Validate simulations against real data distributions
3. **Long-term**: Consider adopting or contributing to community simulation frameworks

The investment in creating a realistic simulator will pay off through credible method comparisons and publishable validation studies.
