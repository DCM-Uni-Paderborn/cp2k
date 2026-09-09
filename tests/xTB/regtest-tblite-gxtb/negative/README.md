# Negative g-xTB capability and input checks

`CH4_gxtb_kp_symmetry_fused_batch_zero.inp` is expected to terminate with a nonzero status and
`IMAGE_BATCH_SIZE must be positive`. CP2K's regression driver treats every nonzero exit status as a
runtime failure, so this fixture is kept outside `TEST_FILES.toml` and is exercised explicitly by
the qualification manifest command.

`H2O_gxtb_minimal_provider_acceleration.inp` is a provider-dependent negative fixture. With a
minimal `lmseidler-integration` library it must terminate with a nonzero status and
`Requested g-xTB acceleration requires an extended save_tblite provider`. It deliberately selects
`MODE MANUAL` and `EXCHANGE_STREAM_BACKEND BATCHED`: an unavailable explicit request must not
silently fall back to the dense implementation. With an extended `tdkuehne-integration` library
this request is supported, so the fixture is not an unconditional regression failure test.
The ordinary `H2O_gxtb_kp_gamma.inp` exercises the minimal provider's default `AUTO` fallback.
