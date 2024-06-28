test_that('Test falciparum switch produces same', {
  parameters_def <- get_parameters()
  parameters_fal <- get_parameters(parasite = "falciparum")
  expect_identical(parameters_def, parameters_fal)
})

test_that('Test vivax model runs', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  sim_res <- run_simulation(timesteps = 100, parameters = vivax_parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test parasite = vivax sets parameters$parasite to vivax', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  expect_identical(vivax_parameters$parasite, "vivax")
})

test_that('Test difference between falciparum and vivax parameter lists', {
  falciparum_parameters <- get_parameters(parasite = "falciparum")
  vivax_parameters <- get_parameters(parasite = "vivax")
  
  expect_true(all(names(falciparum_parameters)[!names(falciparum_parameters) %in% names(vivax_parameters)] %in%
                    c("du","rvm","rva","rb","b0","b1","ib0","kb","theta0","theta1","kv","fv0","av","gammav","iv0","fd0","ad","gammad","d1","id0","kd","ub","uv","gamma1","pvm","init_ivm","init_ib","init_iva")))
  expect_true(all(names(vivax_parameters[!names(vivax_parameters) %in% names(falciparum_parameters)]) %in%
                    c("dpcr_max","dpcr_min","kpcr","apcr50","init_iaa","init_iam","ra","ua","b","philm_min","philm_max","klm","alm50","ca","f","gammal","init_hyp","kmax")))
})

test_that('Test age structure should not change vivax infectivity', {
  falc_parameters <- get_parameters(
    overrides = list(
      human_population = 1,
      init_id  = 0.5))
  
  vivax_parameters <- get_parameters(
    parasite = "vivax",
    overrides = list(
      human_population = 1,
      init_id  = 0.5))
  
  state_mock <- mockery::mock('A', cycle = T)
  mockery::stub(create_variables, 'initial_state', state_mock)
  state_vivax_mock <- mockery::mock(list(human_states = 'A',
                                         hypnozoites_v = 0), cycle = T)
  mockery::stub(create_variables, 'initial_state_vivax', state_vivax_mock)
  
  ages_mock <- mockery::mock(365, cycle = T)
  mockery::stub(create_variables, 'calculate_initial_ages', ages_mock)
  
  falc_variables <- create_variables(falc_parameters)
  vivax_variables <- create_variables(vivax_parameters)
  
  expect_equal(falc_variables$infectivity$get_values(), 0.06761596)
  expect_equal(vivax_variables$infectivity$get_values(), 0.1)
  
  ages_mock <- mockery::mock(365*70, cycle = T)
  mockery::stub(create_variables, 'calculate_initial_ages', ages_mock)
  
  falc_variables <- create_variables(falc_parameters)
  vivax_variables <- create_variables(vivax_parameters)
  
  expect_equal(falc_variables$infectivity$get_values(), 0.03785879)
  expect_equal(vivax_variables$infectivity$get_values(), 0.1)
  
})
