/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#ifndef ACADOS_SOLVER_OCL_Pendulum_H_
#define ACADOS_SOLVER_OCL_Pendulum_H_

#include "acados/utils/types.h"

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#define OCL_PENDULUM_NX     4
#define OCL_PENDULUM_NZ     0
#define OCL_PENDULUM_NU     1
#define OCL_PENDULUM_NP     0
#define OCL_PENDULUM_NBX    1
#define OCL_PENDULUM_NBX0   4
#define OCL_PENDULUM_NBU    1
#define OCL_PENDULUM_NSBX   0
#define OCL_PENDULUM_NSBU   0
#define OCL_PENDULUM_NSH    0
#define OCL_PENDULUM_NSG    0
#define OCL_PENDULUM_NSPHI  0
#define OCL_PENDULUM_NSHN   0
#define OCL_PENDULUM_NSGN   0
#define OCL_PENDULUM_NSPHIN 0
#define OCL_PENDULUM_NSBXN  0
#define OCL_PENDULUM_NS     0
#define OCL_PENDULUM_NSN    0
#define OCL_PENDULUM_NG     0
#define OCL_PENDULUM_NBXN   0
#define OCL_PENDULUM_NGN    0
#define OCL_PENDULUM_NY0    0
#define OCL_PENDULUM_NY     0
#define OCL_PENDULUM_NYN    0
#define OCL_PENDULUM_N      400
#define OCL_PENDULUM_NH     0
#define OCL_PENDULUM_NPHI   0
#define OCL_PENDULUM_NHN    0
#define OCL_PENDULUM_NPHIN  0
#define OCL_PENDULUM_NR     0

#ifdef __cplusplus
extern "C" {
#endif


// ** capsule for solver data **
typedef struct OCL_Pendulum_solver_capsule
{
    // acados objects
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    ocp_nlp_out *sens_out;
    ocp_nlp_solver *nlp_solver;
    void *nlp_opts;
    ocp_nlp_plan_t *nlp_solver_plan;
    ocp_nlp_config *nlp_config;
    ocp_nlp_dims *nlp_dims;

    // number of expected runtime parameters
    unsigned int nlp_np;

    /* external functions */
    // dynamics

    external_function_param_casadi *impl_dae_fun;
    external_function_param_casadi *impl_dae_fun_jac_x_xdot_z;
    external_function_param_casadi *impl_dae_jac_x_xdot_u_z;




    // cost

    external_function_param_casadi *ext_cost_fun;
    external_function_param_casadi *ext_cost_fun_jac;
    external_function_param_casadi *ext_cost_fun_jac_hess;



    external_function_param_casadi ext_cost_0_fun;
    external_function_param_casadi ext_cost_0_fun_jac;
    external_function_param_casadi ext_cost_0_fun_jac_hess;


    external_function_param_casadi ext_cost_e_fun;
    external_function_param_casadi ext_cost_e_fun_jac;
    external_function_param_casadi ext_cost_e_fun_jac_hess;

    // constraints




} OCL_Pendulum_solver_capsule;

ACADOS_SYMBOL_EXPORT OCL_Pendulum_solver_capsule * OCL_Pendulum_acados_create_capsule(void);
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_free_capsule(OCL_Pendulum_solver_capsule *capsule);

ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_create(OCL_Pendulum_solver_capsule * capsule);

ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_reset(OCL_Pendulum_solver_capsule* capsule, int reset_qp_solver_mem);

/**
 * Generic version of OCL_Pendulum_acados_create which allows to use a different number of shooting intervals than
 * the number used for code generation. If new_time_steps=NULL and n_time_steps matches the number used for code
 * generation, the time-steps from code generation is used.
 */
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_create_with_discretization(OCL_Pendulum_solver_capsule * capsule, int n_time_steps, double* new_time_steps);
/**
 * Update the time step vector. Number N must be identical to the currently set number of shooting nodes in the
 * nlp_solver_plan. Returns 0 if no error occurred and a otherwise a value other than 0.
 */
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_update_time_steps(OCL_Pendulum_solver_capsule * capsule, int N, double* new_time_steps);
/**
 * This function is used for updating an already initialized solver with a different number of qp_cond_N.
 */
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_update_qp_solver_cond_N(OCL_Pendulum_solver_capsule * capsule, int qp_solver_cond_N);
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_update_params(OCL_Pendulum_solver_capsule * capsule, int stage, double *value, int np);
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_update_params_sparse(OCL_Pendulum_solver_capsule * capsule, int stage, int *idx, double *p, int n_update);

ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_solve(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_free(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void OCL_Pendulum_acados_print_stats(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int OCL_Pendulum_acados_custom_update(OCL_Pendulum_solver_capsule* capsule, double* data, int data_len);


ACADOS_SYMBOL_EXPORT ocp_nlp_in *OCL_Pendulum_acados_get_nlp_in(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *OCL_Pendulum_acados_get_nlp_out(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *OCL_Pendulum_acados_get_sens_out(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_solver *OCL_Pendulum_acados_get_nlp_solver(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_config *OCL_Pendulum_acados_get_nlp_config(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void *OCL_Pendulum_acados_get_nlp_opts(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_dims *OCL_Pendulum_acados_get_nlp_dims(OCL_Pendulum_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_plan_t *OCL_Pendulum_acados_get_nlp_plan(OCL_Pendulum_solver_capsule * capsule);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SOLVER_OCL_Pendulum_H_
