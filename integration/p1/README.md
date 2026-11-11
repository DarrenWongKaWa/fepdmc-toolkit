# P1 adapter

Source transform of the JJ `update_swap` on the pinned FEP-DMC commit:
installs `linear_da_mod.f90`, adds it to the makefile and inserts the two
delayed-acceptance stages. `apply_p1.py` (with `--dry-run`), `verify_p1.py`;
see [docs/p1.md](../../docs/p1.md) and
[docs/p1-adapter-design.md](../../docs/p1-adapter-design.md).
