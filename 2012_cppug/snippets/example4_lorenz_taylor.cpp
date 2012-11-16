taylor_fixed_order< 25 , 3 > taylor_type stepper;

stepper.do_step(
    fusion::make_vector
    (
        sigma * ( arg2 - arg1 ) ,
        R * arg1 - arg2 - arg1 * arg3 ,
        arg1 * arg2 - b * arg3
    ) , x , t , dt );
