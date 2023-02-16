
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(
  	-Wcast-align
	-Wcast-qual
	-Wctor-dtor-privacy
	-Wdisabled-optimization
	-Wformat=2
	-Winit-self
	-Wlogical-op
	# -Wnoexcept
	# -Wold-style-cast
	-Woverloaded-virtual
	-Wredundant-decls
	# -Wshadow
	-Wsign-promo
	-Wstrict-null-sentinel
	-Wstrict-overflow=5
	-Wswitch-default
	-Wundef
	-Wno-unused
   	-Wmissing-declarations
   	#-Weffc++
	#-Wpadding
)
  # has to be removed at a later stage
  # guess a move ctor and use in code will provide
  # some speed up
  # add_definitions(-Werror)

  add_definitions(
  -Wno-deprecated-copy
    -Wno-error=deprecated-copy
  )
  add_definitions(
   -Wmissing-declarations
  -Wno-error=unused-result
  -Wno-error=sign-compare
  -Wno-error=vla -Wno-error=old-style-cast
  -Wno-error=shadow
  -Wno-error=strict-overflow
  -Wno-error=sign-promo
  -Wno-error=format-nonliteral
  -Wno-error=switch-default
  # mad-ng uses it internally, I modified the header
  -Wno-error=vla-parameter
  # C++ implementation part
  # currently a lot of warnings on gtpsa variant implementation
  -Wno-error=noexcept
  )
endif()
