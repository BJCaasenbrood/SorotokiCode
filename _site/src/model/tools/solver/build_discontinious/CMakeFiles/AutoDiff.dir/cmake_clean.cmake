file(REMOVE_RECURSE
  "libAutoDiff.a"
  "libAutoDiff.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/AutoDiff.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
