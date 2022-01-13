#############################################################
##                                                         ##
## This module is for computing potential energies from    ##
## trajectory files in NAMD                                ##
## Made by Ramon Mendoza, ramendoza@uchicago.edu           ##
## Last update January 12, 2022                            ##
##                                                         ##
#############################################################


if {$print_output} {
 set energy_output_file [open $output_file "w"]
}

proc save_callback {labels values} {
 global saved_labels saved_values
 set saved_labels $labels
 set saved_values $values
}

callback save_callback

proc save_array {} {
 global saved_labels saved_values saved_array
 foreach label $saved_labels value $saved_values {
  set saved_array($label) $value
 }
}

coorfile open dcd $trajectory_file

if {$print_output} {
 puts $energy_output_file [format {%3s%30s} "TS" "POTENTIAL"]
}

set ts_multiplier $ts

while { ![coorfile read] } {
 firstTimestep $ts
 run 0
 
 if {$print_output} {
  save_array
  foreach ts $saved_array(TS) potential $saved_array(POTENTIAL) {
   puts $energy_output_file [format {%4d%#30.*f} $ts $energy_precision $potential] 
  }
 }

 incr ts $ts_multiplier
}

coorfile close
