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
                    c("du","rvm","rva","rb","rid","b0","b1","ib0","kb","theta0","theta1","kv","fv0","av","gammav","iv0","fd0","ad","gammad","d1","id0","kd","ub","uv","ud","gamma1","pvm","init_ivm","init_ib","init_iva","init_id")))
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
      human_population = 1))
  
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

## Test anti-parasite immunity subpatent progression functions
test_that('Test anti-parasite immunity function', {
  vivax_parameters <- get_parameters(parasite = "vivax",
                                     overrides = list(s_proportion = 0,
                                                      d_proportion = 0,
                                                      a_proportion = 0,
                                                      u_proportion = 1,
                                                      t_proportion = 0))
  variables <- create_variables(vivax_parameters)
  index <- variables$state$get_index_of("U")
  dpcr_min <- vivax_parameters$dpcr_min
  dpcr_max <- vivax_parameters$dpcr_max
  apcr50 <- vivax_parameters$apcr50
  kpcr <- vivax_parameters$kpcr
  
  expect_equal(object = anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index)),
    expected = rep(dpcr_max,100))
  
  ## Change initial values of ID, and IDM, check they are the same
  variables$iaa <- individual::DoubleVariable$new(1:100)
  UAA_durations <- anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index))
  
  variables$iaa <- individual::DoubleVariable$new(rep(0,100))
  variables$iam <- individual::DoubleVariable$new(1:100)
  UAM_durations <- anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index))
  
  expect_equal(object = UAA_durations, expected = UAM_durations)
  
  ## Check convergence to min_du at high immunity
  variables$iaa <- individual::DoubleVariable$new(rep(1E7,100))
  variables$iam <- individual::DoubleVariable$new(rep(0,100))
  expect_equal(object = anti_parasite_immunity(
    dpcr_min, dpcr_max, apcr50, kpcr,
    variables$iaa$get_values(index),
    variables$iam$get_values(index)),
    expected = rep(dpcr_min,100),
    tolerance = 1E-2)
})

test_that('that vivax patent prevalence rendering works', {
  
  timestep <- 1
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S', 'Tr'),
    c('U', 'A', 'D', 'D', 'D', 'S')
  )
  birth <- individual::IntegerVariable$new(
    -c(3, 4, 5, 1, 11, 6) * 365
  )
  immunity <- individual::DoubleVariable$new(rep(1, 6))
  vivax_parameters <- get_parameters(parasite = "vivax")
  renderer <- mock_render(1)
  
  process <- create_prevelance_renderer(
    state,
    birth,
    immunity,
    vivax_parameters,
    renderer
  )
  
  process(timestep)
  
  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_detect_lm_730_3650',
    2,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_detect_pcr_730_3650',
    3,
    timestep
  )
  
})

test_that('Test default vivax incidence rendering works', {
  
  timestep <- 0
  year <- 365
  birth <- individual::IntegerVariable$new(
    -c(2, 5, 10, 11) * year
  )
  vivax_parameters <- get_parameters(
    parasite = "vivax")
  
  renderer <- mock_render(1)
  incidence_renderer(
    birth,
    renderer,
    individual::Bitset$new(4)$insert(c(1, 2, 4)),
    'inc_patent_',
    c(0, 2) * year,
    c(5, 10) * year,
    timestep
  )
  
  incidence_probability_renderer(
    birth,
    renderer,
    individual::Bitset$new(4)$insert(seq(4)),
    c(.1, .2, .3, .4),
    'inc_patent_',
    c(0, 2) * year,
    c(5, 10) * year,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_inc_patent_0_1825',
    2,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_inc_patent_730_3650',
    2,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    3,
    'p_inc_patent_0_1825',
    0.3,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    4,
    'p_inc_patent_730_3650',
    .6,
    timestep
  )
})

test_that('relapses are recognised with division between bite infections and relapses', {
  timestep <- 50
  parameters <- get_parameters(parasite = "vivax", overrides = list(human_population = 4))
  
  variables <- list(
    state = individual::CategoricalVariable$new(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    infectivity = individual::DoubleVariable$new(rep(0, 4)),
    recovery_rates = individual::DoubleVariable$new(rep(0, 4)),
    drug = individual::DoubleVariable$new(rep(0, 4)),
    drug_time = individual::DoubleVariable$new(rep(-1, 4)),
    iaa = individual::DoubleVariable$new(rep(0, 4)),
    iam = individual::DoubleVariable$new(rep(0, 4)),
    ica = individual::DoubleVariable$new(rep(0, 4)),
    icm = individual::DoubleVariable$new(rep(0, 4)),
    last_boosted_iaa = individual::DoubleVariable$new(rep(-1, 4)),
    last_boosted_ica = individual::DoubleVariable$new(rep(-1, 4)),
    last_eff_pev_timestep = individual::DoubleVariable$new(rep(-1, 4)),
    pev_profile = individual::IntegerVariable$new(rep(-1, 4)),
    hypnozoites = individual::IntegerVariable$new(c(0, 1, 2, 3))
  )
  
  bernoulli_mock <- mockery::mock(c(2, 4), c(1, 2), 2, cycle = TRUE)
  calc_mock <- mockery::mock(individual::Bitset$new(4)$insert(2))
  mockery::stub(infection_outcome_process, 'bernoulli_multi_p', bernoulli_mock, depth = 2)
  mockery::stub(infection_outcome_process, 'calculate_clinical_infections', calc_mock)
  
  renderer <- mock_render(1)
  infected_humans <- individual::Bitset$new(4)$insert(c(1, 2, 3, 4))
  
  infection_outcome_process(
    timestep = timestep,
    infected_humans,
    variables,
    renderer,
    parameters,
    prob = c(rep(0.5, 4)),
    relative_rate = c(rep(0.5, 4))
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_infections',
    4,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_relapses',
    2,
    timestep
  )
  
})
