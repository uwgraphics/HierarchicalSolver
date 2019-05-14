#!/bin/bash
grep "_mm" $1
perl -pi -e 's|ENABLE_SCALAR_IMPLEMENTATION\((.*?;)\)|\1|' $1
perl -pi -e 's|ENABLE_AVX_IMPLEMENTATION\(.*?;\)||' $1
perl -pi -e 's|ENABLE_SSE_IMPLEMENTATION\(.*?;\)||' $1
perl -pi -e 's|__m128|Tn|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_add_ss\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2 + \3;|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_mul_ss\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2 * \3;|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_sub_ss\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2 - \3;|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_max_ss\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = max(\2, \3);|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_xor_ps\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2 ^ \3;|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_or_ps\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2 \| \3;|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_and_ps\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2 & \3;|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_andnot_ps\(([a-zA-Z_][a-zA-Z0-9_]*),([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2.andnot(\3);|' $1
perl -pi -e 's|([a-zA-Z_][a-zA-Z0-9_]*)=_mm_rsqrt_ss\(([a-zA-Z_][a-zA-Z0-9_]*)\);|\1 = \2.rsqrt();|' $1

grep "_mm" $1