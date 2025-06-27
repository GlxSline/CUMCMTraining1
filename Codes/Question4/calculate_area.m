function area_numb = calculate_area(rho, theta, x_E2, y_E2, r_E1E2, k_E1E2, x_E4, y_E4, r_E3E4, k_E2E3, x_E5, y_E5, k, l)
    x = rho * cos(theta);
    y = rho * sin(theta);
    r_AE2 = sqrt((x - x_E2) ^ 2 + (y - y_E2) ^ 2);
    %* I - II
    theta_1_to_2 = acos((2 * r_E1E2 * r_E1E2 - l ^ 2) / (2 * r_E1E2 * r_E1E2));
    %* II - III
    theta_2_to_3 = acos((2 * r_E3E4 * r_E3E4 - l ^ 2) / (2 * r_E3E4 * r_E3E4));
    %* III - IV
    fun_3_to_4 = @(th) r_E1E2 ^ 2 + (r_E1E2 + k * th) ^ 2 - l ^ 2 - 2 * r_E1E2 * (r_E1E2 + k * th) * cos(th);
    theta_E5 = atan2(y_E5, x_E5);
    options = optimoptions('fsolve', 'Display', 'off');
    theta_3_to_4 = fsolve(fun_3_to_4, theta_E5, options);
    theta_3_to_4 = theta_3_to_4 + 4.5 / k - pi;

    if rho > 4.5

        if abs(rho - k * theta) < 1e-6
            area_numb(1) = 1;
            area_numb(2) = 1;

        else
            area_numb(1) = 4;

            if theta > theta_3_to_4
                area_numb(2) = 4;
            else
                area_numb(2) = 3;
            end

        end

    elseif r_AE2 <= r_E1E2
        area_numb(1) = 2;
        k_E2A = (y - y_E2) / (x - x_E2);
        theta_AE1 = atan((k_E2A - k_E1E2) / (1 + k_E2A * k_E1E2));

        if theta_AE1 > theta_1_to_2
            area_numb(2) = 2;
        else
            area_numb(2) = 1;
        end

    else
        area_numb(1) = 3;
        k_E4A = (y - y_E4) / (x - x_E4);
        theta_AE3 = atan((k_E4A - k_E2E3) / (1 + k_E4A * k_E2E3));

        if theta_AE3 > theta_2_to_3
            area_numb(2) = 3;
        else
            area_numb(2) = 2;
        end

    end
