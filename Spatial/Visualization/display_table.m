function display_table(T)

try
    table((1:size(T,1))',num2str(T{:,1}),T{:,2},T{:,6})
catch
    table((1:size(T,1))',(T{:,1}),T{:,2},T{:,6})
end