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
    max_iterations: 200
    kspace: 150
    output_level: 0
    muelu_xml_file_name: ../muelu_p3.xml
    recompute_preconditioner: no
  # - name: solve_cont
  #   type: hypre
  #   method: hypre_gmres
  #   preconditioner: boomerAMG
  #   tolerance: 1e-5
  #   max_iterations: 200
  #   kspace: 5
  #   output_level: 0

realms:

  - name: realm_1
    mesh: hump2newtop_noplenumZ103x28_2D_p3_ndtw.exo
    use_edges: no
    automatic_decomposition_type: rcb
    polynomial_order: 3
    support_inconsistent_multi_state_restart: yes

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
          velocity: [34.6,0.0]
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
        turbulent_ke: 1e-16
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

        - input_variables_from_file:
            minimum_distance_to_wall: minimum_distance_to_wall

        - turbulence_model_constants:
            SDRWallFactor: 0.625

    turbulence_averaging:
      time_filter_interval: 0.7

      specifications:

        - name: one
          target_name: [Unspecified-2-QUAD]
          compute_reynolds_stress: yes

    data_probes:

      output_frequency: 100

      search_method: stk_octree
      search_tolerance: 1.0e-3
      search_expansion_factor: 2.0

      specifications:
        # - name: probe_bottomwall
        #   from_target_part: bottomwall

        #   line_of_site_specifications:
        #     - name: results_p3/probe_bottomwall
        #       number_of_points: 500
        #       tip_coordinates: [-6.39, 0.0]
        #       tail_coordinates: [4.0, 0.0]

        #   output_variables:
        #     - field_name: tau_wall
        #       field_size: 1
        #     - field_name: pressure
        #       field_size: 1

        - name: probe_profile0
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile0
              number_of_points: 200
              tip_coordinates: [-2.14, 0.0]
              tail_coordinates: [-2.14, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile1
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile1
              number_of_points: 200
              tip_coordinates: [0.65, 0.116101]
              tail_coordinates: [0.65, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile2
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile2
              number_of_points: 200
              tip_coordinates: [0.66, 0.112975]
              tail_coordinates: [0.66, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile3
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile3
              number_of_points: 200
              tip_coordinates: [0.8, 0.0245493]
              tail_coordinates: [0.8, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile4
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile4
              number_of_points: 200
              tip_coordinates: [0.9, 0.00476345]
              tail_coordinates: [0.9, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile5
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile5
              number_of_points: 200
              tip_coordinates: [1.0, 0.0]
              tail_coordinates: [1.0, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile6
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile6
              number_of_points: 200
              tip_coordinates: [1.1, 0.0]
              tail_coordinates: [1.1, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile7
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile7
              number_of_points: 200
              tip_coordinates: [1.2, 0.0]
              tail_coordinates: [1.2, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3

        - name: probe_profile8
          from_target_part: Unspecified-2-QUAD

          line_of_site_specifications:
            - name: results_p3/probe_profile8
              number_of_points: 200
              tip_coordinates: [1.3, 0.0]
              tail_coordinates: [1.3, 0.9]

          output_variables:
            - field_name: velocity
              field_size: 2
            - field_name: reynolds_stress
              field_size: 3


    post_processing:

    - type: surface
      physics: surface_force_and_moment
      output_file_name: results_p3/wallHump.dat
      frequency: 100
      parameters: [0,0]
      target_name: bottomwall

    output:
      output_data_base_name: results_p3/wallHump.e
      output_frequency: 100
      output_node_set: no
      output_variables:
       - velocity
       - pressure
       - pressure_force
       - tau_wall
       - turbulent_ke
       - specific_dissipation_rate
       - minimum_distance_to_wall
       - sst_f_one_blending
       - turbulent_viscosity

    restart:
      restart_data_base_name: restart_p3/wallHump.rst
      restart_frequency: 1000
      restart_time: 0

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      time_step: 1.0e-6
      termination_time: 1.0
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: yes

      realms:
        - realm_1
