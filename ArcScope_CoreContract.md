# ArcScope Core Contract

**Status**: Immutable Architecture Constitution
**Version**: 1.0
**Date**: 2025-12-13

---

## Preamble

This document defines the **non-negotiable architectural constraints** of the ArcScope system. All code, at every layer, must comply with these rules. Violations are not merely bugs—they are **design failures**.

---

## I. Time & Sampling Contract

### R1: Core Time Axis is 1 Hz

**Invariant**: The Core Engine operates on a **per-second time axis**.

- **Sample k** represents the interval `[k, k+1)` seconds
- **Canonical time value**: `t_k = k + 0.5`
- **Film duration**: If a film is `D` seconds long, it has exactly `D` samples (indexed `0` to `D-1`)

### R2: Bridge Output Constraint

**The Bridge layer MUST output per-second raw observations.**

- ✅ **Allowed**: Extract one sample per second from video/audio
- ❌ **Forbidden**: Output sub-second samples, then average in Bridge
- ❌ **Forbidden**: Output pre-fused curves (e.g., `arousalProxy`, `dialogueDensity`)

**Rationale**: Temporal resolution is a global design decision. Changing it affects:
- Database schema
- FilmEngine algorithms (PCA, smoothing, segmentation)
- UI rendering assumptions

### R3: Database Time Invariant

**All `curves` rows MUST have `fps = 1.0`.**

- **Derived constraint**: `length = duration_seconds` (integer)
- **Read operations**: Assume `times[i] = i + 0.5`
- **Write operations**: Never write `fps ≠ 1.0`

**Enforcement**: SQLiteStore SHALL validate this invariant on write.

---

## II. Data Semantic Boundary

### The Three Layers

| Layer | Responsibility | Forbidden Actions |
|-------|---------------|-------------------|
| **Bridge** | Raw observation extraction | Fusion / Diagnosis / Interpretation |
| **FilmEngine** | Fusion / Normalization / Structure | Apple APIs / File I/O / Database access |
| **Store** | Persistence / Query | Inference / Computation |

### R4: Observation vs. Derived Quantities

**A. Raw Observations (Bridge only)**

These are **measurements of the physical world**:

```
✅ motionAmplitude          (pixel displacement)
✅ audioRms                 (audio amplitude)
✅ brightness               (luma average)
✅ saturation               (chroma intensity)
✅ warmth                   (red/blue ratio)
✅ grayness                 (chroma absence)
✅ avgHue                   (dominant color angle)
✅ subtitleDensity          (words/sec)
✅ cutDensity               (cuts/sec)
✅ audioTransient           (onset density)
✅ spectralBalance          (freq ratio)
✅ expositionProbability    (subtitle-derived classifier)
✅ facePresence             (face area ratio)
✅ faceArousal              (CoreML facial AU)
✅ faceValence              (CoreML facial AU)
```

**B. Derived Fusion Quantities (FilmEngine only)**

These are **computed from observations**:

```
❌ infoDensity              (subtitle + visual fusion)
❌ dialogueDensity          (audio + subtitle fusion)
❌ arousalProxy             (motion + audio + face fusion)
❌ pace                     (PCA of motion + cuts + audio)
❌ tension                  (arousal + sound divergence)
```

**C. Semantic / Structural (FilmEngine only)**

These are **interpretations**:

```
❌ SceneType                (classifier output)
❌ Issue                    (diagnostic rule)
❌ ActBoundary              (structure detection)
```

### R5: Enforcement Rules

**Bridge layer violations**:

If `ArcScopeBridgeSample` contains any field from category **B** or **C**, the code is **non-compliant**.

**FilmEngine violations**:

If `FilmEngine.cpp` contains:
- `AVFoundation` imports
- `sqlite3_*` calls
- File I/O (`fopen`, `std::ofstream`)

The code is **non-compliant**.

**Store violations**:

If `SQLiteStore.cpp` contains:
- Mathematical transformations (PCA, smoothing, segmentation)
- Decision logic (if-then rules on data values)

The code is **non-compliant**.

---

## III. ID, Ownership, and Lifecycle

### R6: Film Identity Duality

**Two distinct IDs exist**:

| ID Type | Scope | Type | Stability |
|---------|-------|------|-----------|
| `film_logical_id` | User-facing | String | Permanent (e.g., "blade_runner_1982") |
| `film_pk` | Database-internal | Integer | Volatile (auto-increment) |

**Rules**:

1. **Bridge / UI / URL**: Use `film_logical_id` only
2. **SQLite joins**: Use `film_pk` only
3. **Re-import behavior**: Same `film_logical_id` → overwrite existing `film_pk`

**Current violation**: `films.id` is used as both logical ID and primary key.

**Fix required**:
```sql
ALTER TABLE films ADD COLUMN logical_id TEXT UNIQUE NOT NULL;
CREATE UNIQUE INDEX idx_films_logical_id ON films(logical_id);
```

### R7: Memory Ownership

**Universal Rule**: Whoever allocates memory MUST provide the deallocation interface.

**Bridge layer**:

```c
// ✅ Correct
ArcScopeBridgeCurve* curve = bridge_create_curve(...);
bridge_free_curve(curve);

// ❌ Wrong
free(curve);  // Breaks ABI if Bridge uses new[]
```

**Swift/C++ boundary**:

- **C → Swift**: C layer provides `_free()` function
- **Swift → C**: Swift passes raw pointers, C does NOT own them
- **C++ → C**: C++ allocates with `new[]`, C API uses `delete[]` in `_free()`

**Current state**: ✅ Mostly compliant, but **not documented as a rule**.

### R8: Curve Data Blob Format

**Invariant**: All `curves.data_blob` are stored as **float32 arrays**.

```c
// Write
std::vector<float> values = ...;
sqlite3_bind_blob(stmt, ..., values.data(), values.size() * sizeof(float), ...);

// Read
size_t count = blob_size / sizeof(float);
float* values = (float*)blob_data;
```

**Rationale**:
- 32-bit float is sufficient for UI rendering
- Reduces storage by 50% vs. 64-bit double
- Matches Metal shader expectations

**Current compliance**: ✅ Implemented correctly in SQLiteStore.cpp

---

## IV. Validation & Enforcement

### How to Check Compliance

1. **Time Contract**:
   ```bash
   sqlite3 arcscope.db "SELECT DISTINCT fps FROM curves;"
   # Must output: 1.0 only
   ```

2. **Semantic Boundary**:
   ```bash
   grep -r "arousalProxy\|infoDensity\|dialogueDensity" bridge/include/
   # Must output: nothing
   ```

3. **ID Duality**:
   ```bash
   sqlite3 arcscope.db "PRAGMA table_info(films);"
   # Must have both: id (INTEGER PRIMARY KEY) and logical_id (TEXT UNIQUE)
   ```

4. **Memory Ownership**:
   ```bash
   grep -r "free(" bridge/src/
   # Must output: only bridge_free_* functions
   ```

---

## V. Amendment Process

**This document is immutable** except by unanimous consensus of:
- System architect
- Lead engineer
- Database maintainer

**Rationale for changes** must demonstrate:
1. Current rule is provably incorrect
2. Proposed change does not cascade to other layers
3. Migration path exists for existing data

---

## VI. Penalties for Violations

**During code review**:
- Violation of Time Contract → **Block merge**
- Violation of Semantic Boundary → **Block merge**
- Violation of Memory Ownership → **Block merge**

**In production**:
- Violations create **technical debt** that compounds
- Future changes become exponentially harder
- Data migrations become necessary

**Prevention is mandatory**.

---

## Appendices

### A. Current Non-Compliances

As of 2025-12-13, the following violations exist:

1. ❌ `ArcScopeBridgeSample` contains `info_density`, `dialogue_density`, `arousal_proxy`
   - **Fix**: Remove these fields immediately

2. ❌ `films` table uses `id` as both logical ID and primary key
   - **Fix**: Add `logical_id TEXT UNIQUE` column, migrate data

3. ⚠️ Memory ownership rules are implicit, not documented
   - **Fix**: This document + comments in `arcscope_bridge.h`

### B. Historical Context

**Why 1 Hz?**

- Human perception of film rhythm operates at 0.5–10 second timescales
- Sub-second precision creates noise, not signal
- Storage and computation costs scale linearly with sampling rate

**Why strict layer boundaries?**

- Bridge using FilmEngine logic → Platform lock-in (Apple-only)
- FilmEngine using Store → Violates single responsibility
- Mixing layers → Impossible to unit test

**Why two IDs?**

- Users want stable references ("blade_runner_1982")
- Databases want auto-increment for performance
- Re-imports should replace, not duplicate

---

**End of Core Contract**
