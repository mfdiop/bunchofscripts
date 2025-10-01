# Documentation Index

Complete guide to the IBD Simulation Framework for P. falciparum.

## 📚 Documentation Hierarchy

```
START HERE → README.md
    ↓
    ├─→ New User? → QUICKSTART.md (5 min)
    ├─→ Full Details? → IBD_SIMULATION_README.md
    ├─→ Technical? → ARCHITECTURE.md
    └─→ Status? → STATUS.md
```

## 🚀 Getting Started

### For First-Time Users
1. **[README.md](README.md)** - Project overview and quick start
2. **[QUICKSTART.md](QUICKSTART.md)** - Get running in 5 minutes
3. **[STATUS.md](STATUS.md)** - Current project status

### For Developers
1. **[ARCHITECTURE.md](ARCHITECTURE.md)** - System design and technical details
2. **[IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)** - What was built and how

### For API Reference
1. **[IBD_SIMULATION_README.md](IBD_SIMULATION_README.md)** - Complete documentation

## 📄 Document Descriptions

### README.md
**Purpose**: Main project entry point  
**Audience**: Everyone  
**Contains**:
- Project overview
- Quick start instructions
- Feature list
- Example usage
- Installation guide

**Read this**: To understand what the project does and get started quickly

---

### QUICKSTART.md
**Purpose**: Rapid onboarding  
**Audience**: New users  
**Contains**:
- Installation steps
- First simulation in 3 commands
- Common use cases
- Troubleshooting
- Next steps

**Read this**: To run your first simulation in 5 minutes

---

### IBD_SIMULATION_README.md
**Purpose**: Complete API reference  
**Audience**: Regular users, developers  
**Contains**:
- Full API documentation
- Detailed class/method descriptions
- Advanced usage examples
- P. falciparum parameters
- Extension points
- Scientific background

**Read this**: For comprehensive understanding and advanced usage

---

### ARCHITECTURE.md
**Purpose**: Technical design documentation  
**Audience**: Developers, contributors  
**Contains**:
- System architecture diagrams
- Data flow descriptions
- Algorithm details
- Performance characteristics
- Extension guidelines
- Testing strategy

**Read this**: To understand how the system works internally

---

### IMPLEMENTATION_SUMMARY.md
**Purpose**: What was built and delivered  
**Audience**: Project managers, users  
**Contains**:
- Deliverables checklist
- Feature status
- Code statistics
- Validation results
- Usage examples
- Comparison to alternatives

**Read this**: To understand what was delivered and its capabilities

---

### STATUS.md
**Purpose**: Current project status  
**Audience**: Everyone  
**Contains**:
- Completion checklist
- Quick commands
- File inventory
- Success criteria
- Known limitations

**Read this**: To verify project completion and readiness

---

### INDEX.md (This File)
**Purpose**: Navigation guide  
**Audience**: Everyone  
**Contains**:
- Document hierarchy
- Quick reference guide
- Where to find information

**Read this**: When you need to find specific information

## 🔧 Code Files

### Core Implementation

#### ibd_simulation.py (24 KB)
**Purpose**: Core simulation engine  
**Key Classes**:
- `Segment` - IBD segment with founder label
- `Chromosome` - List of segments
- `Individual` - Diploid individual
- `Population` - Population management
- `ChromosomeStructure` - P. falciparum genome
- `Simulator` - Multi-population simulation

**Use this**: Direct access to simulation primitives

---

#### ibd_calling.py (19 KB)
**Purpose**: IBD detection methods  
**Key Classes**:
- `ThresholdIBDCaller` - Sliding window method
- `HMMIBDCaller` - Hidden Markov Model
- `WindowIBDCaller` - Overlapping windows
- `evaluate_ibd_calling()` - Performance metrics

**Use this**: To implement or use IBD calling methods

---

#### simulation_runner.py (17 KB)
**Purpose**: High-level workflow interface  
**Key Classes**:
- `SimulationRunner` - Complete workflow management
- `run_simple_example()` - Quick demo function

**Use this**: For most simulation tasks (recommended)

---

#### example_simulations.py (12 KB)
**Purpose**: Pre-configured scenarios  
**Key Functions**:
- `example_1_basic_simulation()` - Single population
- `example_2_with_bottleneck()` - Bottleneck event
- `example_3_two_populations_with_migration()` - Migration
- `example_4_isolated_populations()` - No migration
- `example_5_varying_snp_density()` - SNP comparison
- `run_all_examples()` - Batch execution

**Use this**: As templates for your simulations

---

#### test_basic_functionality.py (9 KB)
**Purpose**: Unit tests  
**Key Functions**:
- Tests for all core components
- Validation checks
- `run_all_tests()` - Execute all tests

**Use this**: To verify installation and functionality

---

#### requirements.txt
**Purpose**: Python dependencies  
**Contains**: numpy, pandas, scipy, matplotlib, seaborn

**Use this**: With `pip install -r requirements.txt`

## 🎯 Quick Reference: Where to Find...

### Installation Instructions
- **README.md** - Quick install
- **QUICKSTART.md** - Detailed install

### First Example
- **QUICKSTART.md** - Simplest example
- **README.md** - Basic usage
- **simulation_runner.py** - `run_simple_example()`

### API Reference
- **IBD_SIMULATION_README.md** - Complete API
- Code docstrings - Inline documentation

### How It Works
- **ARCHITECTURE.md** - Technical details
- **IMPLEMENTATION_SUMMARY.md** - High-level overview

### Examples
- **example_simulations.py** - 5 scenarios
- **IBD_SIMULATION_README.md** - Usage examples
- **QUICKSTART.md** - Common use cases

### Extending the Framework
- **ARCHITECTURE.md** - Extension points
- **IBD_SIMULATION_README.md** - Advanced usage
- Code comments - Implementation details

### Performance
- **ARCHITECTURE.md** - Performance characteristics
- **IMPLEMENTATION_SUMMARY.md** - Tested configurations
- **STATUS.md** - Performance notes

### Validation
- **test_basic_functionality.py** - Unit tests
- **IMPLEMENTATION_SUMMARY.md** - Validation results
- **STATUS.md** - Validation checklist

### Troubleshooting
- **QUICKSTART.md** - Common issues
- **STATUS.md** - Known limitations

### Scientific Background
- **README.md** - Overview
- **IBD_SIMULATION_README.md** - Detailed background
- **ARCHITECTURE.md** - Algorithm details

## 📖 Reading Paths

### Path 1: Quick Start (15 minutes)
1. README.md - Overview (3 min)
2. QUICKSTART.md - Get started (5 min)
3. Run `test_basic_functionality.py` (2 min)
4. Run simple example (5 min)

### Path 2: Regular User (1 hour)
1. README.md - Overview (5 min)
2. QUICKSTART.md - Examples (10 min)
3. IBD_SIMULATION_README.md - API (30 min)
4. example_simulations.py - Code review (15 min)

### Path 3: Developer (2-3 hours)
1. README.md - Overview (5 min)
2. IMPLEMENTATION_SUMMARY.md - What was built (15 min)
3. ARCHITECTURE.md - How it works (45 min)
4. Code files - Implementation review (60-90 min)
5. test_basic_functionality.py - Tests (15 min)

### Path 4: Research Use (ongoing)
1. QUICKSTART.md - Get started (10 min)
2. example_simulations.py - Run examples (30 min)
3. Modify for your needs (ongoing)
4. IBD_SIMULATION_README.md - Reference as needed

## 🔍 Search Guide

### I want to...

**Install the framework**
→ README.md or QUICKSTART.md

**Run my first simulation**
→ QUICKSTART.md

**Understand all features**
→ README.md or IMPLEMENTATION_SUMMARY.md

**Learn the API**
→ IBD_SIMULATION_README.md

**Understand the algorithms**
→ ARCHITECTURE.md

**See example code**
→ example_simulations.py or QUICKSTART.md

**Add a new IBD method**
→ ARCHITECTURE.md (Extension Points) + ibd_calling.py

**Add demographic events**
→ IBD_SIMULATION_README.md (Advanced Usage) + ibd_simulation.py

**Compare to other tools**
→ IMPLEMENTATION_SUMMARY.md (Comparison section)

**Optimize performance**
→ ARCHITECTURE.md (Performance section) + STATUS.md

**Verify it works**
→ test_basic_functionality.py + STATUS.md

**Understand IBD in malaria**
→ README.md or IBD_SIMULATION_README.md (Scientific Background)

## 📊 Document Statistics

| Document | Size | Lines | Purpose |
|----------|------|-------|---------|
| README.md | 6 KB | 190 | Main overview |
| QUICKSTART.md | 5 KB | 210 | Quick start |
| IBD_SIMULATION_README.md | 12 KB | 430 | Full API |
| ARCHITECTURE.md | 13 KB | 500 | Technical design |
| IMPLEMENTATION_SUMMARY.md | 10 KB | 350 | Deliverables |
| STATUS.md | 4 KB | 150 | Current status |
| INDEX.md | 6 KB | 250 | This file |

**Total Documentation**: ~56 KB, ~2,080 lines

## 🎓 Learning Sequence

### Beginner Level
1. Read README.md
2. Follow QUICKSTART.md
3. Run examples
4. Modify parameters

### Intermediate Level
1. Read IBD_SIMULATION_README.md
2. Study example_simulations.py
3. Create custom scenarios
4. Implement simple extensions

### Advanced Level
1. Read ARCHITECTURE.md
2. Study core code files
3. Add new IBD methods
4. Optimize performance
5. Contribute improvements

## ✅ Completion Status

- [x] All code files complete
- [x] All documentation complete
- [x] All tests passing
- [x] Examples working
- [x] Installation verified

**Project Status**: ✅ COMPLETE AND READY

---

**Last Updated**: October 1, 2025  
**Maintained By**: Project maintainers  
**Purpose**: Navigate the IBD Simulation Framework documentation
