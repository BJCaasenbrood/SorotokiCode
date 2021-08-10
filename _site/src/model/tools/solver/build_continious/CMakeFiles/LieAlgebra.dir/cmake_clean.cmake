file(REMOVE_RECURSE
  "libLieAlgebra.a"
  "libLieAlgebra.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/LieAlgebra.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
