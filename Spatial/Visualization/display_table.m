function display_table(T)

table((1:size(T,1))',num2str(T{:,1}),T{:,2},T{:,6})