test_suite ising

test ggg
 assert_real_equal(0.0,viscosity(0.0))

end test

test gg
 integer:: state(5,5)
 real::energy
 state=1
 call get_energy(state,energy)
 assert_real_equal(2.0,energy)
end test

end test_suite
