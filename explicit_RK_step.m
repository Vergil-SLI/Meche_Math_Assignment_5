function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    k = zeros(length(XA), length(BT_struct.B));
    sum = 0;
    num_evals = 0;
    %h_list = logspace(h(1), h(2), h(3));
    %h_avg = (tspan(2) - tspan(1))

    for i = 1:length(BT_struct.B)
        
        t_input = t + BT_struct.C(i) * h;
        X_input = 0;       
        for j = 1:i-1
            X_input = X_input + BT_struct.A(i, j) * k(:,j );
        end
        X_input = XA + h*X_input;

        k(:, i) = rate_func_in(t_input, X_input);

        num_evals = num_evals + 1;

        sum = sum + BT_struct.B(1,i) * k(:, i);
    end
 
    XB = XA + h*sum;
end