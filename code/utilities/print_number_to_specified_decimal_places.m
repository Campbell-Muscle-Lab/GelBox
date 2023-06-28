function text_string=print_number_to_specifified_decimal_places( ...
    input_number,decimal_places)

% Create the format string
format_string=sprintf('%%0.%if',decimal_places);

% Now create the string to be evaluated
eval_string=sprintf('text_string=sprintf(''%s'',%f);', ...
    format_string,input_number);

% And evaluate
eval(eval_string)

