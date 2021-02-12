function [input1a, input1b] = get_flows_for_pool(flow_in,flow_out,pool, show)
% note: pool 0-3
input1a = flow_in(:,pool+1);
input1b = flow_out(:,pool);

if show
    plot(input1a)
    hold on;
    plot(input1b)
    title("Water flows");
    legend("q_{in}", "q_{out}")
end

end

