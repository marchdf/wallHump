# -*- mode: yaml -*-
Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0
    muelu_xml_file_name: ../../muelu.xml

realms:

  - name: realm_1
    mesh: hump2newtop_noplenumZ205x55_2D.exo
    use_edges: no
    automatic_decomposition_type: rcb
    polynomial_order: 2

    time_step_control:
     target_courant: 10.0
     time_step_change_factor: 1.2

    equation_systems:
      name: theEqSys
      max_iterations: 3

      solver_system_specification:
        velocity: solve_scalar
        turbulent_ke: solve_scalar
        specific_dissipation_rate: solve_scalar
        pressure: solve_cont

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - ShearStressTransport:
            name: mySST
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:
      - constant: ic_1
        target_name: [Unspecified-2-QUAD]
        value:
          pressure: 0
          velocity: [34.6,0.0,0.0]
          turbulent_ke: 0.00108
          specific_dissipation_rate: 7710.9

    material_properties:
      target_name: [Unspecified-2-QUAD]
      specifications:
        - name: density
          type: constant
          value: 1.185
        - name: viscosity
          type: constant
          value: 1.8398e-5

    boundary_conditions:

    - wall_boundary_condition: bc_wall
      target_name: bottomwall
      wall_user_data:
        velocity: [0,0]
        turbulent_ke: 0.0
        use_wall_function: no

    - symmetry_boundary_condition: bc_symTop
      target_name: top
      symmetry_user_data:

    - inflow_boundary_condition: bc_inflow
      target_name: inlet
      inflow_user_data:
        velocity: [34.6,0.0]
        turbulent_ke: 0.00108
        specific_dissipation_rate: 7710.9

    - open_boundary_condition: bc_open
      target_name: outlet
      open_user_data:
        velocity: [0,0]
        pressure: 0.0
        turbulent_ke: 0.00108
        specific_dissipation_rate: 7710.9

    solution_options:
      name: myOptions
      turbulence_model: sst

      options:
        - hybrid_factor:
            velocity: 0.0
            turbulent_ke: 1.0
            specific_dissipation_rate: 1.0

        - alpha:
            velocity: 0.0

        - limiter:
            pressure: no
            velocity: no
            turbulent_ke: no
            specific_dissipation_rate: yes

        - projected_nodal_gradient:
            velocity: element
            pressure: element
            turbulent_ke: element
            specific_dissipation_rate: element

        - turbulence_model_constants:
            SDRWallFactor: 0.625

    restart:
      restart_data_base_name: hump2newtop_noplenumZ205x55_2D_p2.exo
      restart_start: 0
      restart_frequency: 1

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 0
      time_step: 1.0e-6
      time_stepping_type: fixed 
      time_step_count: 0
      second_order_accuracy: no

      realms:
        - realm_1
