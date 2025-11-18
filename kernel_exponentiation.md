# Kernel exponentiation

## Function implementing exponentiation
- `kernel/propagator.m` (`P=propagator(spin_system,L,timestep)`) computes the propagator `exp((-1i)*L*timestep)` via a scaled-and-squared Taylor series (see lines 1-239) with optional GPU acceleration and caching.

## How the exponentiation is carried out
- Small matrices shortcut: if `size(L,1)<spin_system.tols.small_matrix`, it returns `expm((-1i*timestep)*L)` directly (kernel/propagator.m:36-39).
- Propagator caching: when `'prop_cache'` is enabled, a hash of `{L,timestep,prop_chop}` drives load/save of cached `P` in `spin_system.sys.scratch` (kernel/propagator.m:41-91, 243-272).
- Generator setup: builds `A=clean_up(spin_system,(-1i*timestep)*L,prop_chop)` and reports density/sparsity before exponentiation (kernel/propagator.m:93-100).
- Norm guardrails: uses `cheap_norm(A)`; if `>1e9` it errors; if `>1024` it tightens `prop_chop` to `eps(''double'')` and warns; if `>16` it warns (kernel/propagator.m:102-123).
- Scaling selection: chooses `n_squarings=max([0 ceil(log2(mat_norm))])`, scales `A` by `1/2^n` and plans to square `P` back `n` times (kernel/propagator.m:125-136).
- Taylor series core (exponentiation):
  - GPU path when `'gpu'` enabled and matrix size >500: converts `A` to `gpuArray`, initialises `P=speye` and `next_term=gpuArray.speye`; loop forms each term until `nnz(next_term)==0`, multiplying as `(1/n)*A*next_term` for sparse `A` or `(1/n)*next_term*A` otherwise, cleaning small elements each step (kernel/propagator.m:140-164). After summation it gathers only the accumulated `next_term` additions to CPU and reports iterations (kernel/propagator.m:165-181).
  - CPU path otherwise: same Taylor loop on CPU with `speye`, switching multiplication order based on `issparse(A)` (kernel/propagator.m:159-182).
- Post-series clean-up: clears temporaries, cleans `P` with `clean_up`, reports density (kernel/propagator.m:185-195).
- Squaring stage: if `n_squarings>0`, repeatedly squares `P` (`P=clean_up(P*P,prop_chop)`) either on GPU (with gather at end) or CPU, reporting each step and final density (kernel/propagator.m:197-239).
- Optional cache save: writes the propagator to disk when caching is enabled and the computation took >0.01 s, using MAT v7.3 (kernel/propagator.m:243-272).

All exponentiation in the kernel is handled by `kernel/propagator.m`; other propagation helpers call this routine to form the exponential propagator when needed.
