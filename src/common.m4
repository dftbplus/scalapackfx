dnl
dnl Undefining some M4 builtins to avoid conflicts with Fortran code
dnl invoke them via the builtin() command if needed.
dnl
undefine(`len')dnl
undefine(`index')dnl
undefine(`shift')dnl


dnl Sets a variable ($1) to the value of an optional argument ($2)
dnl if present or to a default value ($3) otherwise.
dnl
define(`_handle_inoptflag',`dnl
if (present($2)) then
  $1 = $2
else
  $1 = $3
end if
')


dnl Sets an optional output argument ($1) if present to a certain value ($2).
dnl
define(`_handle_outoptflag', `dnl
if (present($1)) then
  $1 = $2
end if
')


dnl Moves allocation from an array ($1) to an optional target ($2).
dnl
define(`_optmovealloc', `dnl
if (present($2)) then
  call move_alloc($1, $2)
end if
')


dnl Allocates an array ($1) to a minimal size ($2) with an actual size
dnl stored in $3. If the optional allocatable argument ($4) is present
dnl and big enough, its allocation transfer will be transfered instead of
dnl a new allocation.
dnl
define(`_move_minoptalloc', `dnl
if (present($4)) then
  if (size($4) >= $2) then
    call move_alloc($4, $1)
  else
    deallocate($4)
  end if
end if
if (.not. allocated($1)) then
  allocate($1($2))
end if
$3 = size($1)
')


