program template

  intrinsic              :: selected_real_kind
  integer,  parameter    :: wp = selected_real_kind(15)
  real(wp), allocatable  :: identity(:,:)
  integer                :: i, n_dim_identity

  n_dim_identity = 2
    identity(1:n_dim_identity, 1:n_dim_identity) = 0.0d0
    forall(i = 1:n_dim_identity) identity(i, i) = 1.0d0

end program template
