from affine_gaps import needleman_wunsch_gotoh_alignment

insulin = "GIVEQCCTSICSLYQLENYCN"
glucagon = "HSQGTFTSDYSKYLDSRAEQDFV"
aligned_insulin, aligned_glucagon, aligned_score = needleman_wunsch_gotoh_alignment(insulin, glucagon)

print("Alignment 1:", aligned_insulin)  # GI-V---EQCC-TSICSLY---QL-ENYCN-
print("Alignment 2:", aligned_glucagon) # --D-FVHSQGTFTSDYSKYLDSRAEQDF--V
print("Score:", aligned_score)          # 41