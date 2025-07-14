function tandemWings_Tester
    figure;
    names = {};
    hold on;
    
    for x = 0.00:0.001:0.006
        % disp(x)
        try
            [x_path, y_path] = tandemWings_Simluate(x);
        catch
            disp(["Skipped", num2str(x)])
        end
        plot(x_path, y_path, "LineWidth",1);
        names{end + 1} = "Nose mass: " + num2str(x);
    end
    hold off;
    xlabel("X-Position (m)");
    ylabel("Y-Position (m)");
    title("Tandem Wing Flight with Varying C.G.'s");
    grid on;
    
    legend(names);
end