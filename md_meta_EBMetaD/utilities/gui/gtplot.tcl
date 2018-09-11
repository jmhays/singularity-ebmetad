#
# GT plot v 1.0
# 
# Gareth Tribello
# 

package provide gtPlot 1.0

namespace eval ::gtPlot:: {
    namespace export gtPlot
    variable maxgrid  100 ; # Maximum number of points in a grid that you allow to plot (depends on how long you are willing to wait)
    variable ntics_d    3 ; # Default number of major tics
    variable appearance   ; # Controls how the graph looks
    variable bounds       ; # The automatic boundaries for the graph 
    variable zoomID       ; # The device used by zoom
    variable fesID        ; # Stuff used by fes
    variable scrollID     ; # Stuff used by scrolling tool
    variable pwindow      ; # The window in which we plot 
    variable plot         ; # This tells the graph when to plot
    variable kT           ; # The value of kT
    variable axisd        ; # A handle to the axis properties dialogue window 
    variable canv         ; # A handle to the window in which we draw the graph
}

# This gets arguments from a list
proc gtPlot::getargs {arglist tag defl {n 1}} {
   set pos [lsearch $arglist $tag]
   if {$pos<0}  { return $defl}
   return [join [ lindex $arglist [expr $pos+$n ] ] ]
}    

proc gtPlot::setDefaults {args} {
   variable appearance
   variable kT

   # These control the appearance of the axis
   set appearance(xfixtics) 0
   set appearance(yfixtics) 0
   set appearance(xmajticl) 4.0
   set appearance(ymajticl) 4.0
   set appearance(xminticl) 2.0
   set appearance(yminticl) 2.0
   set appearance(xnmintic) 1
   set appearance(ynmintic) 1
   set appearance(zmajticl) 4.0
   set appearance(zminticl) 2.0
   set appearance(znmintic) 1
   # Set defaults for kT (300K and kJ mol-1)
   set kT(T)                300
   set kT(units)            0.00831447
   set kT(kT)               [ expr $kT(T)*$kT(units) ] 
}

# This creates the graph
proc ::gtPlot::create {wind args} {
    variable appearance
    variable fesID
    variable bounds
    variable plot
    variable canv
    
    # Interpret the arguments
    set size [getargs $args "-size" {}]

    # Setup some default stuff
    setDefaults
    
    # Store the appearance information
    set appearance(font) [getargs $args "-font" {}]
    set appearance(pad) [getargs $args "-padding" {}]
#    set appearance(trajcolors) [getargs $args "-trajcol" {}]
#    set appearance(ncolors) [getargs $args "-nfescolors" {}]
#    set fescolors [getargs $args "-fescol" {}]    

    # Creates a canvas with a border
    frame $wind
    set pad $appearance(pad)
    eval { canvas $wind.canvas -highlightthickness 0 -borderwidth 0 -background white } \
      -width [ expr $size + 2 * $pad ] -height [ expr $size + 2 * $pad ]
   
    # Store the location of the canvas
    set canv $wind.canvas
 
    # Setup the grid stuff
    grid $wind.canvas -sticky news
    grid rowconfigure $wind 0 -weight 1
    grid columnconfigure $wind 0 -weight 1

    # This will set up our zoom tool
#    bind $wind.canvas <Double-Button-1> [namespace code { gtPlot::autoAxes %W} ]             ; # Return to automatic axis
    bind $wind.canvas <Shift-Button-1> [namespace code { gtPlot::ZoomStart %W %x %y} ]       ; # Set the initial point
    bind $wind.canvas <Shift-B1-Motion> [namespace code { gtPlot::ZoomMove %W %x %y} ]       ; # Adjust the box size 
    bind $wind.canvas <Shift-B1-ButtonRelease> [namespace code { gtPlot::ZoomEnd %W %x %y} ] ; # Delete the zoom box
    # This is for the free energy calculation tool
    set fesID(st) 0                                                                             ; # This is so free energy differences work 
    bind $wind.canvas <Control-Button-1> [namespace code { gtPlot::FesStart %W %x %y} ]         ; # Create a fes integration box
    bind $wind.canvas <Control-B1-Motion> [namespace code { gtPlot::FesMove %W %x %y} ]         ; # Adjust the box size
    bind $wind.canvas <Control-B1-ButtonRelease> [namespace code { gtPlot::FesEnd %W %x %y} ]   ; # Calculate the free energy 
    bind $wind.canvas <Button-1> [namespace code { gtPlot::FesClear %W %x %y} ]                 ; # Clear all the fes data
    # This is for the little thing that keeps track of our position
    bind $wind.canvas <Motion> [namespace code { gtPlot::writeCoord %W %x %y} ] 
    # This changes the size of the graph if the window changes size 
    bind $wind.canvas <Configure> [namespace code { gtPlot::resizeWindow } ] 

    # These traces all control the tic appearance
    trace add variable appearance(xmajticl) write [namespace code "drawTics -axis x"] 
    trace add variable appearance(ymajticl) write [namespace code "drawTics -axis y"] 
    trace add variable appearance(xminticl) write [namespace code "drawTics -axis x"]
    trace add variable appearance(yminticl) write [namespace code "drawTics -axis y"] 
    trace add variable appearance(xnmintic) write [namespace code "drawTics -axis x"] 
    trace add variable appearance(ynmintic) write [namespace code "drawTics -axis y"] 
  
    return 
}

# This destroys the graph
proc ::gtPlot::destroy {args} {
    variable appearance
    trace remove variable appearance(xmajticl) write [namespace code "drawTics -axis x"]      
    trace remove variable appearance(ymajticl) write [namespace code "drawTics -axis y"]      
    trace remove variable appearance(xminticl) write [namespace code "drawTics -axis x"]
    trace remove variable appearance(yminticl) write [namespace code "drawTics -axis y"]
    trace remove variable appearance(xnmintic) write [namespace code "drawTics -axis x"]         
    trace remove variable appearance(ynmintic) write [namespace code "drawTics -axis y"]   
}

# This removes data from the plot
proc ::gtPlot::remove {args} {
    variable data 
    variable canv

    set tag [getargs $args "-datatag" {}]     ;# What data are we deleting

    # Check there is data plotted with this tag
    if { [llength [$canv find withtag $tag] ]==0 } { return }
 
    $canv delete $tag
    $canv delete basinsize
    
    # Keep track of what data is displayed
    set data([ list $tag displayed ]) 0 

    # Remove the integrated free energy data if this is a fes
    if { $data([list $tag isfes])==1 } { $canv delete fesint } 

    # Check if any data is left (if there is none delete the axis)
    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
            set tag $data($keyitem)
            # Check if data is displayed
            if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
       }
    }
    if { $flag == 0 } { 
       $canv delete axis 
       $canv delete whereIam
    }
}

proc gtPlot::axisDialogueDestroy {args} {
  variable bounds
  trace remove variable bounds(xticspacing) write [namespace code "drawTics -axis x"]
  trace remove variable bounds(yticspacing) write [namespace code "drawTics -axis y"]
}

proc ::gtPlot::axisPropertiesDialogue {args} {
  variable axisd
  variable appearance
  variable bounds
  variable data

  # Check that there is data plotted
  set flag 0
  foreach keyitem [array names data] {
     # Look for all the tag names
     if { [ lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          # Check if data is displayed
          if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
     }
  }     
  if { $flag == 0 } { return }

  if [winfo exists .axis] {
    wm deiconify $axisd
    return
  }

  set bounds(txmin) $bounds(xmin)  ; set bounds(tymin) $bounds(ymin)
  set bounds(txmax) $bounds(xmax)  ; set bounds(tymax) $bounds(ymax)

  set axisd [toplevel ".axis"]
  wm title $axisd "Axis Properties"
  wm resizable $axisd 0 0
  bind $axisd <Destroy> gtPlot::axisDialogueDestroy

  frame $axisd.top -padx 1m -pady 1m
  frame $axisd.cont -padx 1m -pady 1m

  labelframe $axisd.x -relief ridge -bd 2 -text "x axis" -padx 2m -pady 2m
  label $axisd.x.minl -text "Min" -anchor e
  entry $axisd.x.min -width 5 -textvariable gtPlot::bounds(txmin)
  label $axisd.x.maxl -text "Max" -anchor e
  entry $axisd.x.max -width 5 -textvariable gtPlot::bounds(txmax)
  grid $axisd.x.minl -row 1 -column 1 -sticky e
  grid $axisd.x.min -row 1 -column 2 -sticky e
  grid $axisd.x.maxl -row 1 -column 3 -sticky e
  grid $axisd.x.max -row 1 -column 4 -sticky e 

  labelframe $axisd.x.majfram -relief ridge -bd 1 -text "Major Tics" -padx 1m -pady 1m
  label $axisd.x.majfram.spal -text "spacing" -anchor e
  entry $axisd.x.majfram.spa -width 5 -textvariable gtPlot::bounds(xticspacing)
  label $axisd.x.majfram.lenl -text "tic length" -anchor e
  eval spinbox $axisd.x.majfram.len -textvariable gtPlot::appearance(xmajticl) -from 1.0 -to 8.0 -increment 0.25 -width 5
  label $axisd.x.majfram.fixtl -text "Fix tic spacing" -anchor e
  checkbutton $axisd.x.majfram.fixt -variable gtPlot::appearance(xfixtics)
  grid $axisd.x.majfram.spal -row 1 -column 1 -sticky e
  grid $axisd.x.majfram.spa  -row 1 -column 2 -sticky e
  grid $axisd.x.majfram.lenl -row 1 -column 3 -sticky e
  grid $axisd.x.majfram.len  -row 1 -column 4 -sticky e
  grid $axisd.x.majfram.fixtl -row 2 -column 2 -columnspan 2 -sticky e
  grid $axisd.x.majfram.fixt -row 2 -column 4 -sticky e

  labelframe $axisd.x.minfram -relief ridge -bd 1 -text "Minor Tics" -padx 1m -pady 1m
  label $axisd.x.minfram.nl -text "n. tics" -anchor e
  eval spinbox $axisd.x.minfram.n -textvariable gtPlot::appearance(xnmintic) -from 0 -to 5 -increment 1 -width 5
  label $axisd.x.minfram.lenl -text "tic length" -anchor e
  eval spinbox $axisd.x.minfram.len -textvariable gtPlot::appearance(xminticl) -from 1.0 -to 4.0 -increment 0.125 -width 5
  grid $axisd.x.minfram.nl -row 1 -column 1 -sticky e
  grid $axisd.x.minfram.n -row 1 -column 2 -sticky e
  grid $axisd.x.minfram.lenl -row 1 -column 3 -sticky e
  grid $axisd.x.minfram.len -row 1 -column 4 -sticky e

  grid $axisd.x.majfram -row 2 -column 1 -columnspan 4 -sticky e
  grid $axisd.x.minfram -row 3 -column 1 -columnspan 4 -sticky e
  pack $axisd.x -in $axisd.cont -side left -fill both

  labelframe $axisd.y -relief ridge -bd 2 -text "y axis" -padx 2m -pady 2m
  label $axisd.y.minl -text "Min" -anchor e
  entry $axisd.y.min -width 5 -textvariable gtPlot::bounds(tymin)
  label $axisd.y.maxl -text "Max" -anchor e
  entry $axisd.y.max -width 5 -textvariable gtPlot::bounds(tymax)
  grid $axisd.y.minl -row 1 -column 1 -sticky e
  grid $axisd.y.min -row 1 -column 2 -sticky e
  grid $axisd.y.maxl -row 1 -column 3 -sticky e
  grid $axisd.y.max -row 1 -column 4 -sticky e

  labelframe $axisd.y.majfram -relief ridge -bd 1 -text "Major Tics" -padx 1m -pady 1m
  label $axisd.y.majfram.spal -text "spacing" -anchor e
  entry $axisd.y.majfram.spa -width 5 -textvariable gtPlot::bounds(yticspacing)
  label $axisd.y.majfram.lenl -text "tic length" -anchor e
  eval spinbox $axisd.y.majfram.len -textvariable gtPlot::appearance(ymajticl) -from 1.0 -to 8.0 -increment 0.25 -width 5
  label $axisd.y.majfram.fixtl -text "Fix tic spacing" -anchor e
  checkbutton $axisd.y.majfram.fixt -variable gtPlot::appearance(yfixtics)
  grid $axisd.y.majfram.spal -row 1 -column 1 -sticky e
  grid $axisd.y.majfram.spa  -row 1 -column 2 -sticky e
  grid $axisd.y.majfram.lenl -row 1 -column 3 -sticky e
  grid $axisd.y.majfram.len  -row 1 -column 4 -sticky e
  grid $axisd.y.majfram.fixtl -row 2 -column 2 -columnspan 2 -sticky e
  grid $axisd.y.majfram.fixt -row 2 -column 4 -sticky e
  
  labelframe $axisd.y.minfram -relief ridge -bd 1 -text "Minor Tics" -padx 1m -pady 1m
  label $axisd.y.minfram.nl -text "n. tics" -anchor e
  eval spinbox $axisd.y.minfram.n -textvariable gtPlot::appearance(ynmintic) -from 0 -to 5 -increment 1 -width 5
  label $axisd.y.minfram.lenl -text "tic length" -anchor e
  eval spinbox $axisd.y.minfram.len -textvariable gtPlot::appearance(yminticl) -from 1.0 -to 4.0 -increment 0.125 -width 5
  grid $axisd.y.minfram.nl -row 1 -column 1 -sticky e
  grid $axisd.y.minfram.n -row 1 -column 2 -sticky e
  grid $axisd.y.minfram.lenl -row 1 -column 3 -sticky e
  grid $axisd.y.minfram.len -row 1 -column 4 -sticky e

  grid $axisd.y.majfram -row 2 -column 1 -columnspan 4 -sticky e
  grid $axisd.y.minfram -row 3 -column 1 -columnspan 4 -sticky e
  pack $axisd.y -in $axisd.cont -side right -fill both

  pack $axisd.cont -in $axisd.top -side top -fill both
  frame $axisd.but -padx 3m -pady 3m
  pack [ button $axisd.but.r  -text "Replot" -relief raised -command { [namespace code gtPlot::setAxis]; [namespace code gtPlot::resizeWindow] } ] -in $axisd.but -side right
  pack [ button $axisd.but.d  -text "Dismiss" -relief raised -command { [namespace code gtPlot::axisDialogueDestroy]; destroy .axis} ] -in $axisd.but -side right
  pack $axisd.but -in $axisd.top -side bottom -fill both
  pack $axisd.top -fill both

  trace add variable bounds(xticspacing) write [namespace code "drawTics -axis x"]
  trace add variable bounds(yticspacing) write [namespace code "drawTics -axis y"] 
}

proc ::gtPlot::setkT {args} {
    variable kT

    set tmp [getargs $args "-kT" 0]

    if { $tmp==0 } {
       set kT(T) 300
       set kT(units) 0

       set t [toplevel ".kt"]
       wm title $t "Set value of kT"
       wm resizable $t 0 0

       frame $t.top -padx 2m -pady 2m

       frame $t.up
       label $t.up.tlab -text "Temperature" -anchor e
       entry $t.up.t -width 5 -textvariable gtPlot::kT(T)
       label $t.up.ulab -text "Units" -anchor e
       menubutton $t.up.u -relief raised -bd 2 -direction flush -width 5 \
                -textvariable gtPlot::kT(units) -menu $t.up.u.menu
       menu $t.up.u.menu
       $t.up.u.menu delete 0 end
       $t.up.u configure -state disabled

       $t.up.u configure -state normal
       $t.up.u.menu add radiobutton -label "kJ/mol" -value 0.00831447 -variable gtPlot::kT(units)
       $t.up.u configure -state normal
       $t.up.u.menu add radiobutton -label "kcal/mol" -value 0.001987191 -variable gtPlot::kT(units)
       $t.up.u configure -state normal
       $t.up.u.menu add radiobutton -label "eV" -value 0.00008617343 -variable gtPlot::kT(units)
       $t.up.u configure -state normal
       $t.up.u.menu add radiobutton -label "dlpoly" -value 0.831451115 -variable gtPlot::kT(units)
       $t.up.u configure -state normal
       $t.up.u.menu add radiobutton -label "Rydbergs" -value 0.0000063363125 -variable gtPlot::kT(units)   
       $t.up.u configure -state normal
       $t.up.u.menu add radiobutton -label "Natural" -value 1.0 -variable gtPlot::kT(units)
 
       grid $t.up.tlab -row 1 -column 1 -sticky e
       grid $t.up.t -row 1 -column 2 -sticky e
       grid $t.up.ulab -row 1 -column 3 -sticky e
       grid $t.up.u -row 1 -column 4 -sticky e 

       pack $t.up -in $t.top -side top -fill both
       pack [ button $t.ok -text "OK" -relief raised -command { set gtPlot::kT(kT) [expr $gtPlot::kT(T)*$gtPlot::kT(units)]; destroy .kt } ] -in $t.top -side bottom
       pack $t.top -fill both
    } else {
       set kT(kT) $tmp
    }
}

proc ::gtPlot::readData {args} {
    variable data
    variable appearance

    # Get the number of data points
    set tag [getargs $args "-datatag" {}]
    set data([list $tag np]) [getargs $args "-npoints" {}]
    set data([list $tag isfes]) [getargs $args "-isfes" 0]

    set check [lsearch $args "-nonewbounds"]
    set data([list $tag fes fixedz]) 0

    set nsets 0   ; set safe 0
    # Store the tag names 
    foreach keyitem [array names data] {
       # The != ensures we overwrite any data with the same tag
       if { [ lsearch $keyitem tag ] != -1 && $data($keyitem) != $tag } { incr nsets }
       # If we are using old bounds for z then we check that there was old data
       if { $check>0 &&  $data($keyitem)==$tag } { set safe 1 }
    }
    set data([ list tag $nsets ]) $tag

    # use zbounds from an old calculation  
    if { $safe==1 } { set data([list $tag fes fixedz]) 1 }

    # Get x data
    set pos [lsearch $args "-xdata"]
    if {$pos<0} { tk_messageBox -icon error -type ok -title Message -message "No x data";  return }   
    set data([list $tag x]) [ lindex $args [expr $pos+1] ]
    getRange -tag $tag -axis x

    # Get y data 
    set pos [lsearch $args "-ydata"]
    if {$pos<0} { tk_messageBox -icon error -type ok -title Message -message "No y data";  return }   
    set data([list $tag y]) [ lindex $args [expr $pos+1] ]
    getRange -tag $tag -axis y

    # Get z data (if its there)
    set pos [lsearch $args "-zdata"]
    if {$pos>0} { 
       set data([list $tag is3d]) 1
       set data([list $tag z]) [ lindex $args [expr $pos+1] ] 
       # Get zrange if we are not using one from previous calculation
       if { $data([list $tag fes fixedz]) == 0 } { getRange -tag $tag -axis z }
       # Fix the zscale (i.e. don't adjust with zooming) 
       set pos2 [lsearch $args "-fixedz"]
       if { $pos2>0 } { set data([list $tag fes fixedz]) 1 }
       getRange -tag $tag -axis z
       # Put the read in data onto a grid
       setupGrid -tag $tag
    } else {
       set data([list $tag is3d]) 0
    }

    # Data is currently not displayed
    set data([list $tag displayed]) 0
    # Currently there are no highlights
    set data([list $tag arehighlights]) 0
    # We want to see every point
    set data([ list $tag stride ]) 1
}

proc gtPlot::getRange {args} {
    variable data
  
    set axis [getargs $args "-axis" {}]
    set tag [getargs $args "-tag" {}]
    set np $data([list $tag np])   ; # Retrieve number of data points

    # We have to check the data has the correct number of elements
    if { [llength $data([ list $tag $axis])]!=$np } {
       tk_messageBox -icon error -type ok -title Message -message "Not enough $axis data"
       return
    }

    for { set i 0} { $i<$np } { incr i } {
        set p [lindex $data([list $tag $axis]) $i]
	# Find maximum and minimum in the data
        if { $i == 0 } { 
	    set min $p;  set max $p;
	} 
	if { $p > $max } { set max $p }
	if { $p < $min } { set min $p } 
    }
    
    if { $min == $max } { 
        # If the data is all on one point reset the x and y boundaries artificially
        set data([ list $tag $axis min ]) [ expr $min - 1 ]; set data([ list $tag $axis max ]) [ expr $max + 1 ]; 
    } else {
        # Otherwise change the boundaries by 5% of the range
        set data([ list $tag $axis min ]) [ expr $min - ( 0.05 *( $max - $min ) ) ]  ;
        set data([ list $tag $axis max ]) [ expr $max + ( 0.05 *( $max - $min ) ) ]  ;  
    }
}

proc gtPlot::setupGrid {args} {
    variable data

    set tag [getargs $args "-tag" {}]

    # Find the number of points in the grid in the x direction
    set kk 0
    while { [lindex $data([list $tag x]) $kk] == [lindex $data([list $tag x]) [expr $kk+1] ] } { incr kk }
    incr kk

    # Setup the grid
    set np $data([list $tag np])
    set i 0
    for { set nn 0 } { $nn<$np } { set nn [expr $nn + $kk] } {
        for { set j 0 } { $j<$kk } { incr j } {
            set data([list $tag fes $i $j]) \
                [list [lindex $data([list $tag x]) [expr $nn+$j] ] \
                      [lindex $data([list $tag y]) [expr $nn+$j] ] \
                      [lindex $data([list $tag z]) [expr $nn+$j] ] ]
        }
        incr i
    }

    # We have now transfered the data onto a grid so we 
    # can get rid of the data we read in
    unset data([list $tag x])
    unset data([list $tag y])
    unset data([list $tag z])

    # Setup x and y grid size
    set data([list $tag y ngrid]) $kk
    set data([list $tag x ngrid]) $i
    # Setup x and y pixel widths
    set data([list $tag x pixwidth]) [expr [lindex $data([list $tag fes 1 0]) 0] - [lindex $data([list $tag fes 0 0]) 0] ]
    set data([list $tag y pixwidth]) [expr [lindex $data([list $tag fes 0 1]) 1] - [lindex $data([list $tag fes 0 0]) 1] ]

    # Some debugging - lets test how we do with the grid
    #for { set nn 0 } { $nn<$i } { incr nn } {
    #    for { set j 0 } { $j<$kk } { incr j } {
    #         puts "Grid point at gx, gy $nn $j x [lindex $data([list $tag fes $nn $j]) 0] y [lindex $data([list $tag fes $nn $j]) 1] z [lindex $data([list $tag fes $nn $j]) 2]"
    #    }
    #}
}

proc ::gtPlot::setHighlights {args} {
    variable appearance
    variable data

    # This allows us to highlight certain data points
    set tag [getargs $args "-datatag" {}]
    set data([list $tag arehighlights]) 1
    set data([list $tag highlight center]) [getargs $args "-highlight" {}]
    set data([list $tag highlight isline]) [getargs $args "-lines" {}]
    set data([list $tag highlight nbefore]) [getargs $args "-nbefore" {}]
    set data([list $tag highlight nafter]) [getargs $args "-nafter" {}]
    set appearance(trajcolors) [getargs $args "-colors" {}]
}

proc ::gtPlot::updateHighlight {args} {
    variable appearance
    variable data
    variable canv

    set tag [getargs $args "-datatag" {}]
    set data([list $tag highlight center]) [getargs $args "-highlight" {}]
    set data([list $tag highlight isline]) [getargs $args "-lines" {}]
    set data([list $tag highlight nbefore]) [getargs $args "-nbefore" {}]
    set data([list $tag highlight nafter]) [getargs $args "-nafter" {}]
    set appearance(trajcolors) [getargs $args "-colors" {}]

    # Get rid of highlighting of drawn points
    foreach id [$canv find withtag highlighted] {
       # Recover the frame number this corresponds to
       set taglist [$canv gettags $id]
       set listindex [lsearch -glob $taglist $tag*]
       set frametag [lindex $taglist $listindex]
       lassign [split $frametag :] foo num

       # Check if point is plotted (if so reconfigure and dtag
       if { [expr $num%$data([list $tag stride])] == 0 } {
          $canv itemconfigure $id -outline black -fill white
          $canv dtag $id highlighted
       } 
    }

    # Get rid of old highlight stuff
    $canv delete highlighted
    $canv delete highlight_line

    drawHighlights $canv -datatag $tag
}

proc gtPlot::drawHighlights {args} {
    variable appearance
    variable bounds
    variable data
    variable canv

    set tag [getargs $args "-datatag" {}]

    # Setup a colorscale
    if { $data([list $tag highlight nbefore])==0 && $data([list $tag highlight nafter])==0 } {
       set colormap [createColormap -colorscale $appearance(trajcolors) -ncolors 2]
    } elseif { $data([list $tag highlight nbefore])>$data([list $tag highlight nafter ]) } {
       set colormap [createColormap -colorscale $appearance(trajcolors) -ncolors [expr $data([list $tag highlight nbefore]) + 1] ]
    } else {
       set colormap [createColormap -colorscale $appearance(trajcolors) -ncolors [expr $data([list $tag highlight nafter]) + 1]]
    } 

    # Add highlight on central point
    if { [expr $data([list $tag highlight center])%$data([list $tag stride])] == 0 } {
       set id [$canv find withtag $tag$data([list $tag highlight center]) ]
       $canv itemconfigure $id -outline black -fill [lindex $colormap 0]
       $canv raise $id displayed
       $canv addtag highlighted withtag $tag$data([list $tag highlight center]) 
    } else {
       set x [lindex $data([list $tag x]) $data([list $tag highlight center])]
       set y [lindex $data([list $tag y]) $data([list $tag highlight center])]
       if { $x >= $bounds(xmin) && $x <= $bounds(xmax) && $y >= $bounds(ymin) && $y <= $bounds(ymax) } {
           drawPixel -at $x $y -tag "pixel notfes $tag$data([list $tag highlight center]) $tag highlighted" -fill [lindex $colormap 0] -outline black
       } else {
           drawPixel -at $x $y -tag "pixel notfes $tag$data([list $tag highlight center]) $tag highlighted" -fill [lindex $colormap 0]  -outline black
       }
    }

    # Now add highlights on points before and after center  
    for { set i 1 } { $i<=$data([list $tag highlight nbefore]) } { incr i } {
        set frno [expr $data([list $tag highlight center]) - $i ]
        if { $frno>=0 } {
           if { [expr $frno%$data([list $tag stride])] == 0 } {
              set id [$canv find withtag $tag$frno]
              $canv itemconfigure $id -outline black -fill [lindex $colormap $i]
              $canv raise $id displayed
              $canv addtag highlighted withtag $tag$frno
           } else {
              set x [lindex $data([list $tag x]) $frno] 
              set y [lindex $data([list $tag y]) $frno]
              if { $x >= $bounds(xmin) && $x <= $bounds(xmax) && $y >= $bounds(ymin) && $y <= $bounds(ymax) } {
                  drawPixel -at $x $y -tag "pixel notfes $tag$frno $tag highlighted" -fill [lindex $colormap $i] -outline black
              } else {
                  drawPixel -at $x $y -tag "pixel notfes $tag$frno $tag highlighted" -fill [lindex $colormap $i] -outline black
              }
           }

           if { $data([list $tag highlight isline]) == 1 } {
              set grno [expr $frno + 1 ]
              set x1 [lindex $data([list $tag x]) $frno]
              set y1 [lindex $data([list $tag y]) $frno]
              set x2 [lindex $data([list $tag x]) $grno]
              set y2 [lindex $data([list $tag y]) $grno]
              drawLine -from $x1 $y1 -to $x2 $y2 -tag "$tag notfes highlight_line hl:$i" -fill [lindex $colormap $i] -width 0.1
              set lid [$canv find withtag hl:$i]
              $canv lower $lid displayed
           } 
        }
    } 

    for { set i 1 } { $i<=$data([list $tag highlight nafter]) } { incr i } {
        set frno [expr $data([list $tag highlight center]) + $i ]
        if { $frno<$data([list $tag np]) } {
           if { [expr $frno%$data([list $tag stride])] == 0 } {
              set id [$canv find withtag $tag$frno]
              $canv itemconfigure $id -outline black -fill [lindex $colormap $i]
              $canv raise $id displayed
              $canv addtag highlighted withtag $tag$frno
           } else {
              set x [lindex $data([list $tag x]) $frno]
              set y [lindex $data([list $tag y]) $frno]
              if { $x >= $bounds(xmin) && $x <= $bounds(xmax) && $y >= $bounds(ymin) && $y <= $bounds(ymax) } {
                  drawPixel -at $x $y -tag "pixel notfes $tag$frno $tag highlighted" -fill [lindex $colormap $i] -outline black
              } else {
                  drawPixel -at $x $y -tag "pixel notfes $tag$frno $tag highlighted" -fill [lindex $colormap $i] -outline black
              }
           }

           if { $data([list $tag highlight isline]) == 1 } {
              set grno [expr $frno - 1 ]
              set x1 [lindex $data([list $tag x]) $frno]
              set y1 [lindex $data([list $tag y]) $frno]
              set x2 [lindex $data([list $tag x]) $grno]
              set y2 [lindex $data([list $tag y]) $grno]
              drawLine -from $x1 $y1 -to $x2 $y2 -tag "$tag notfes highlight_line al:$i" -fill [lindex $colormap $i] -width 0.1
              set lid [$canv find withtag al:$i]
              $canv lower $lid displayed
           }
        } 
    }  
}

proc ::gtPlot::drawBasinSize {args} {
    variable canv

    set center [getargs $args "-center" {}]
    set shape [getargs $args "-shape" {}]
    set size [getargs $args "-size" {}]

    # Remove any old basin size
    $canv delete basinsize

    # Extract position of center
    set xc [lindex $center 0]  ;   set yc [lindex $center 1]
    # Extract shape   
    set xx [lindex $shape 0]   ;   set xy [lindex $shape 1]
    set yx [lindex $shape 2]   ;   set yy [lindex $shape 3] 

#    puts "Shape:"
#    puts "$xx $xy"
#    puts "$yx $yy"

    if ($yx!=$yx) {
       puts "ERROR - mismatch between off diagonal elements of the covariance"
       return
    }

    # Diagonalize the input 2x2 matrix
    set delta [ expr sqrt( ($xx-$yy)*($xx-$yy) + 4*$xy*$yx ) ]
    set eval1 [ expr ( ($xx+$yy) + $delta ) / 2.0 ]
    set eval2 [ expr ( ($xx+$yy) - $delta ) / 2.0 ]

    # Get the eigenvectors using trick from Aitkins Molecular Quantum mechanics (p34)
    # (assume here that $xy=$yx)
    set theta [ expr -0.5*atan ( 2*$xy / ($yy-$xx) ) ]
    set evec1(0) [ expr cos($theta) ]  ;  set evec1(1) [ expr sin($theta) ]
    set evec2(0) [ expr -sin($theta) ] ;  set evec2(1) [ expr cos($theta) ]

    set pi 3.141592653589793
    set nincrements 50

    set increment [expr 2*$pi / $nincrements]

    set ncircles [llength $size]

    # These are the zero positions of the circles
    for {set j 0} { $j<$ncircles} { incr j } {
       set tvec1 [expr [lindex $size $j]*( sqrt($eval1)*$evec1(0) + sqrt($eval2)*$evec1(1) ) ] 
       set tvec2 [expr [lindex $size $j]*( sqrt($eval1)*$evec2(0) + sqrt($eval2)*$evec2(1) ) ]
       set oldpx($j) [expr $xc + cos(0)*$tvec1 ]
       set oldpy($j) [expr $yc + sin(0)*$tvec2 ] 
    }

    for {set i 1} {$i<=$nincrements} {incr i} {
        set theta [ expr $i*$increment ]
        for {set j 0} { $j<$ncircles} { incr j } {
           set tvec1 [expr [lindex $size $j]*( sqrt($eval1)*$evec1(0) + sqrt($eval2)*$evec1(1) ) ]
           set tvec2 [expr [lindex $size $j]*( sqrt($eval1)*$evec2(0) + sqrt($eval2)*$evec2(1) ) ]
           set px [expr $xc + cos($theta)*$tvec1 ]
           set py [expr $yc + sin($theta)*$tvec2 ]
        
           # Draw a line that will make up one portion of our elipse
           drawLine -from $oldpx($j) $oldpy($j) -to $px $py -tag "notfes basinsize" -fill blue -width 0.1
           # Save the old coordinates
           set oldpx($j) $px    ;  set oldpy($j) $py
        }
    }
}

proc gtPlot::getAutoBound {args} {
    variable data

    set searchtype [getargs $args "-type" {}]
    set axis [getargs $args "-axis" {}]

    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)   
          # Check the data is displayed 
          if { $data([ list $tag displayed ]) == 1 } {
              if { $flag == 0 } {
                   set result $data([ list $tag $axis $searchtype ]) ; incr flag
              } 
              if { $searchtype == "min" && $data([ list $tag $axis $searchtype ])<$result } {
                   set result $data([ list $tag $axis $searchtype ])
              } elseif { $searchtype == "max" && $data([ list $tag $axis $searchtype ])>$result } {
                   set result $data([ list $tag $axis $searchtype ])
              } 
          }
       }
    }
    return $result
}

proc gtPlot::setAxis {args} {
    variable appearance
    variable bounds
    variable ntics_d
    
    # Read axis data
    set xmin [getargs $args "-xrange" {none} 1]
    set xmax [getargs $args "-xrange" {none} 2]
    set ymin [getargs $args "-yrange" {none} 1]
    set ymax [getargs $args "-yrange" {none} 2]

    if { $xmin=="none" || $xmax=="none" || $ymin=="none" || $ymax=="none" } {
       set xmin $bounds(txmin);   set xmax $bounds(txmax)
       set ymin $bounds(tymin);   set ymax $bounds(tymax)
    } 
    
    # Check the axis are sane
    if { $xmax<$xmin || $ymax<$ymin } { return }

    # Now set the axis for real
    set bounds(xmin) $xmin ; set bounds(xmax) $xmax 
    set bounds(ymin) $ymin ; set bounds(ymax) $ymax 

    # Setup other stuff for drawing points from range
    set bounds(xcenter) [expr ( $xmin + $xmax ) / 2.0 ]
    set bounds(ycenter) [expr ( $ymin + $ymax ) / 2.0 ]
    set bounds(xrange) [expr $xmax - $xmin ]
    set bounds(yrange) [expr $ymax - $ymin ]

    # Find the origin (on the canvas)
    findOrigin -axis "x"   ;   findOrigin -axis "y"

    if { $appearance(xfixtics)==0 } {
       if { [expr $xmax - $bounds(xorigin_points) ] > [ expr $bounds(xorigin_points) - $xmin ] } {
          set bounds(xticspacing) [ format %1.0e [ expr ( $xmax - $bounds(xorigin_points) ) / $ntics_d ] ]
       } else {
          set bounds(xticspacing) [ format %1.0e [ expr ( $bounds(xorigin_points) - $xmin ) / $ntics_d ] ]
       } 
    }
    if { $appearance(yfixtics)==0 } {
       if { [expr $ymax - $bounds(yorigin_points) ] > [ expr $bounds(yorigin_points) - $ymin ] } {
          set bounds(yticspacing) [ format %1.0e [ expr ( $ymax - $bounds(yorigin_points) ) / $ntics_d ] ]
       } else {
          set bounds(yticspacing) [ format %1.0e [ expr ( $bounds(yorigin_points) - $ymin ) / $ntics_d ] ]
       }
    }
}

proc gtPlot::findOrigin {args} {
    variable appearance
    variable bounds
    variable canv
    
    set axis [getargs $args "-axis" {}]

    if { $axis == "x" } {
	set csize [expr [ winfo width $canv ] - 2.0*$appearance(pad) ]
        set min $bounds(xmin)
        set max $bounds(xmax)
    } else {
	set csize [expr [ winfo height $canv ] - 2.0*$appearance(pad) ]
        set min $bounds(ymin)
        set max $bounds(ymax)
    }
    
    if { $min<0 && $max>0 && $axis=="x" } {
        set bounds(xorigin_points) 0
        set bounds(xorigin_canvas) [ point2canvas $canv -point 0.0 -axis $axis ]
    } elseif { $min<0 && $max>0 && $axis=="y" } {
        set bounds(yorigin_points) 0
        set bounds(yorigin_canvas) [ point2canvas $canv -point 0.0 -axis $axis ]
    } elseif { $max<=0 && $axis=="x" } {
        set bounds(xorigin_points) $max
        set bounds(xorigin_canvas) [ expr $appearance(pad) + $csize ]
    } elseif { $min>=0 && $axis=="x" } {
        set bounds(xorigin_points) $min
        set bounds(xorigin_canvas) $appearance(pad) 
    } elseif { $max<=0 && $axis=="y" } {
        set bounds(yorigin_points) $max
        set bounds(yorigin_canvas) $appearance(pad) 
    } elseif { $min>=0 && $axis=="y" } {
        set bounds(yorigin_points) $min
        set bounds(yorigin_canvas) [ expr $appearance(pad) + $csize ]
    } else {
        tk_messageBox -icon error -type ok -title Message -message "Couldn't work out find origin"    
    }    
}

proc gtPlot::drawAxis {args} {
    variable appearance
    variable bounds
    variable canv
    
    set axis [getargs $args "-axis" {}]
    
    if { $axis == "x" } {
	set csize [expr [ winfo width $canv ] - 2.0*$appearance(pad) ]
	$canv create line $appearance(pad) $bounds(yorigin_canvas) [ expr $csize + $appearance(pad) ] $bounds(yorigin_canvas) -tag "notfes axis"
    } else {
	set csize [expr [ winfo height $canv ] - 2.0*$appearance(pad) ]
	$canv create line $bounds(xorigin_canvas) $appearance(pad) $bounds(xorigin_canvas) [ expr $csize + $appearance(pad) ] -tag "notfes axis"
    }
    drawTics -axis $axis    
}

proc ::gtPlot::getFesScaleAppearance {args} {
   variable appearance 
   return [list $appearance(zticspacing) $appearance(znmintic) $appearance(zmajticl) $appearance(zminticl)]
}

proc ::gtPlot::setFesScaleAppearance {ticspacing nmintics majticl minticl} {
   variable appearance
   set appearance(zticspacing) $ticspacing
   set appearance(znmintic) $nmintics
   set appearance(zmajticl) $majticl
   set appearance(zminticl) $minticl
   createScalez
}

proc gtPlot::createScalez {args} {
   variable appearance
   variable bounds
   variable canv

   set ncolors [getargs $args "-ncolors" {}]
   set fescolors [getargs $args "-colorscale" {}]

   # Check that ticspacing is a double
   if {![string is double -strict $appearance(zticspacing)] } { return }

   $canv delete colorscale

   set colormap [ createColormap -colorscale $fescolors -ncolors $ncolors ]
   set iheight [expr 0.5*[winfo height $canv] / $ncolors ]
   set sstart [expr 0.75*[winfo height $canv] ]

   # Set the x position of the colorscale
   set xl [expr [winfo width $canv] - 0.8*$appearance(pad) ] 
   set xu [expr [winfo width $canv] - 0.6*$appearance(pad) ]

   # Draw the colorscale
   for {set i 0 } { $i<$ncolors } { incr i } {
      set yl [expr $sstart - $i*$iheight] ; set yu [expr $yl - $iheight]
      set color [ lindex $colormap $i ]
      $canv create rectangle $xl $yl $xu $yu -tag "colorscale notfes" -fill $color -outline $color
   }

   # Create a rectangle around the color scale
   set yu [expr $sstart - $ncolors*$iheight ] ;   
   $canv create rectangle $xl $sstart $xu $yu -tag "colorscale notfes" -outline black

   # Now create tics
   set ticspacing $appearance(zticspacing)
   set nmtics $appearance(znmintic)
   set majticl $appearance(zmajticl)
   set minticl $appearance(zminticl)
   set oaxis [expr [winfo width $canv] - 0.6*$appearance(pad) ] 
   set ntics [ expr int( ( $bounds(zmax) - $bounds(zmin) ) / $ticspacing) ]
   set realspacing [expr ( 0.5*[winfo height $canv] / ( $bounds(zmax) - $bounds(zmin) ) ) * $ticspacing ]

   # Draw tics
   for { set i 0 } { $i<=$ntics } { incr i } {
      set lab [ format %.4g [ expr $bounds(zmin) + ( $i * $ticspacing ) ] ]
      set p [ expr $sstart - $i*$realspacing  ]
      $canv create line $oaxis $p [ expr $oaxis + $majticl ] $p -tag "colorscale notfes" 
      $canv create text [ expr $oaxis + 2 * $majticl ] $p -font $appearance(font) -anchor w -text $lab -tag colorscale 

      # Draw minor tics 
      for { set j 1 } { $j <=$nmtics } { incr j } {
          set mp [ expr $p + ($j * $realspacing/($nmtics+1) ) ] 
          $canv create line $oaxis $mp [ expr $oaxis + $minticl ] $mp -tag "colorscale notfes"
      }
   }

}

proc gtPlot::drawTics {args} {
    variable appearance
    variable data
    variable bounds
    variable canv

    # Check that data is plotted
    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
            set tag $data($keyitem)
            # Check if data is displayed
            if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
       }
    }
    if { $flag == 0 } { return }
    
    set axis [getargs $args "-axis" {}]
    $canv delete tics$axis

    if { $axis == "x" } {
       set csize [expr [ winfo width $canv ] - 2.0*$appearance(pad) ]
       set min $bounds(xmin)
       set max $bounds(xmax)
       set ticdir 1
       set zero $bounds(xorigin_points)
       set oaxis $bounds(yorigin_canvas)
       set majticl $appearance(xmajticl)
       set minticl $appearance(xminticl)
       set nmtics $appearance(xnmintic)
       set ticspacing $bounds(xticspacing)
    } else {
       set csize [expr [ winfo height $canv ] - 2.0*$appearance(pad) ]
       set min $bounds(ymin)
       set max $bounds(ymax)
       set ticdir -1
       set zero $bounds(yorigin_points)
       set oaxis $bounds(xorigin_canvas)
       set majticl $appearance(ymajticl)
       set minticl $appearance(yminticl)
       set nmtics $appearance(ynmintic)
       set ticspacing $bounds(yticspacing)
    }

    # Check that major ticspacing is a number
    if { ![string is double -strict $ticspacing] } { return }

    # Draw minor tics between origin and first major tic
    for { set j 1 } { $j <=$nmtics } { incr j } {
       set mlab [ expr $zero + ( $j * $ticspacing/($nmtics+1) ) ]
       set mp [ point2canvas -point $mlab -axis $axis ]
       if { $mlab<$max && $axis=="x" } {
          $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "notfes axis tics$axis"
       } elseif { $mlab<$max && $axis=="y" } {
          $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "notfes axis tics$axis"
       }
    }

    # Draw tics on +ve axis
    set ntics [ expr int( ( $max - $zero ) / $ticspacing) ]
    for { set i 1 } { $i <= $ntics } { incr i } {
       set lab [ format %.4g [ expr $zero + ( $i * $ticspacing ) ] ] 
       set p [ point2canvas -point $lab -axis $axis ]
       if { $lab<$max && $axis=="x" } { 
          $canv create line $p $oaxis $p [ expr $oaxis + $majticl ] -tag "notfes axis tics$axis"
          $canv create text $p [ expr $oaxis + 2 * $majticl ] -font $appearance(font) -anchor n -text $lab -tag "notfes axis tics$axis label"
       } elseif { $lab<$max && $axis=="y" } {
          $canv create line $oaxis $p [ expr $oaxis - $majticl ] $p -tag "notfes axis tics$axis"
          $canv create text [ expr $oaxis - 2 * $majticl ] $p -font $appearance(font) -anchor e -text $lab -tag "notfes axis tics$axis label"
       }
       for { set j 1 } { $j <=$nmtics } { incr j } {
          set mlab [ expr $lab + ( $j * $ticspacing/($nmtics+1) ) ]
          set mp [ point2canvas -point $mlab -axis $axis ]
          if { $mlab<$max && $axis=="x" } { 
             $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "notfes axis tics$axis"
          } elseif { $mlab<$max && $axis=="y" } {
             $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "notfes axis tics$axis"
          }
       }
    }

    # Draw minor tics between origin and first major tic
    for { set j 1 } { $j <=$nmtics } { incr j } {
       set mlab [ expr $zero - ( $j * $ticspacing/($nmtics+1) ) ]
       set mp [ point2canvas -point $mlab -axis $axis ]
       if { $mlab<$max && $axis=="x" } {
          $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "notfes axis tics$axis"
       } elseif { $mlab<$max && $axis=="y" } {
          $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "notfes axis tics$axis"
       }  
    }

    # Draw tics on -ve axis
    set ntics [expr int( ( $zero - $min ) / $ticspacing) ]
    for { set i 1 } { $i <= $ntics } { incr i } {
       set lab [ format %.4g [ expr $zero - ( $i * $ticspacing ) ] ] 
       set p [ point2canvas -point $lab -axis $axis ]
       if { $lab>$min && $axis=="x" } {
          $canv create line $p $oaxis $p [ expr $oaxis + $majticl ] -tag "notfes axis tics$axis"
          $canv create text $p [ expr $oaxis + 2 * $majticl ] -font $appearance(font) -anchor n -text $lab -tag "notfes axis tics$axis label"
       } elseif { $lab>$min && $axis=="y" } {
          $canv create line $oaxis $p [ expr $oaxis - $majticl ] $p -tag "notfes axis tics$axis"
          $canv create text [ expr $oaxis - 2 * $majticl ] $p -font $appearance(font) -anchor e -text $lab -tag "notfes axis tics$axis label"
       }
       for { set j 1 } { $j <=$nmtics } { incr j } {
          set mlab [ expr $lab - ( $j * $ticspacing/($nmtics+1) ) ]
          set mp [ point2canvas -point $mlab -axis $axis ]
          if { $lab>$min && $axis=="x" } {
             $canv create line $mp $oaxis $mp [ expr $oaxis + $minticl ] -tag "notfes axis tics$axis"
          } elseif { $lab>$min && $axis=="y" } {
             $canv create line $oaxis $mp [ expr $oaxis - $minticl ] $mp -tag "notfes axis tics$axis"
          }
       }
    }
}
    
proc ::gtPlot::plot2d {args} {
    variable appearance
    variable data 
    variable bounds
    variable canv   
 
    # Get arguments
    set appearance(size)  [getargs $args "-pointsize" {}]
    #set appearance(fescolors) [getargs $args "-fescolors" {}]
    #set appearance(ncolors) [getargs $args "-ncolors" {}]
    set tag               [getargs $args "-datatag" {}]
    set data([ list $tag fescolors]) [getargs $args "-fescolors" {}]
    set data([ list $tag ncolors]) [getargs $args "-ncolors" {}]
    set data([ list $tag stride ]) [getargs $args "-stride" {}]
    set data([ list $tag color ]) [getargs $args "-color" white]

    # Delete any old plot with this tag
    $canv delete $tag

    # Setup how to plot the data
    set data([list $tag drawpoints]) 0 ; set data([list $tag drawline]) 0 ; set data([list $tag drawbars]) 0
    if { [getargs $args "-drawpoints" 0] != 0 } { set data([list $tag drawpoints]) 1 }
    if { [getargs $args "-drawline" 0] != 0 } { set data([list $tag drawline]) 1 }
    if { [getargs $args "-drawbars" 0] != 0 } { 
         set data([list $tag drawbars]) 1 
         # For bar charts we have to extend the x axis from the range in the data
         set data([ list $tag x max]) [expr $data([ list $tag x max ]) + [lindex $data([list $tag x]) 1] - [lindex $data([list $tag x]) 0] ] 
    }

    # This keeps track of which data is plotted
    set data([ list $tag displayed ]) 1

    # Plot the data with automatic axis
    autoAxes 

    # This updates the free energy difference information if it is there
    if { $data([ list $tag isfes ])==1 } { writeFesData $canv }
}

# This allows us to control the behaviour of the z scale on zooming
proc ::gtPlot::setproperty {args} {
   variable data
   set tag               [getargs $args "-datatag" {}]
   set data([ list $tag fes fixedz ]) [getargs $args "-fesZoomStyle" 0]
}

proc gtPlot::resizeWindow {args} {
    variable data
    variable canv
    
    # This checks if there is any data to plot 
    # if there isn't then we return
    set flag 0
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
       }
    }
    if {$flag==0 } { return }
 
    $canv delete whereIam        ; # Delete the where to plot (its now in the wrong place)
    $canv delete axis            ; # Delete the old axis
    $canv delete pixel           ; # Delete the old data
    $canv delete highlight_line  ; # Delete the old trajectory lines
    $canv delete basinsize       ; # Delete the sizes of any old basins
    $canv delete dataline

    # Find the origin (the one on the canvas changes with resizes)
    findOrigin -axis "x"   ;   findOrigin -axis "y"

    # Draw the new axis here
    drawAxis -axis "x"
    drawAxis -axis "y"

    # This plots all the data sets
    set flag 0
    foreach keyitem [array names data] {
       if { [lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          if { $data([ list $tag displayed ]) == 1 } {
               if { $data([ list $tag drawbars ])==1 } {
                   barplotter -tag $tag   ; # This plots a bar chart
               } elseif { $data([ list $tag is3d ]) == 0 } {
                   2dplotter -tag $tag    ; # This plots data points 
               } elseif { $data([ list $tag is3d ]) == 1 && $flag == 0 } {
                   3dplotter -tag $tag    ; # This plots free energy surfaces  
                   set flag 1
               } else {
                   tk_messageBox -icon error -type ok -title Message -message "Can't plot more than one fes at once"    
                   return
               }
          }
       }
    }
}

proc ::gtPlot::autoAxes {args} {
    variable data
    variable canv

    # Check that there is data to plot
    set flag 0   
    foreach keyitem [array names data] {
       # Look for all the tag names
       if { [ lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
       }
    }
    if { $flag==0 } { return }

    # Delete old axes and old data points
    $canv delete axis
    $canv delete pixel
    $canv delete whereIam
    $canv delete highlight_line
    $canv delete basinsize 
    $canv delete dataline

    # This tells the plotter to set auto axis 
    setAxis -xrange [getAutoBound -type min -axis x] [getAutoBound -type max -axis x] \
            -yrange [getAutoBound -type min -axis y] [getAutoBound -type max -axis y]      

    # Draw the new axis here
    drawAxis -axis "x"
    drawAxis -axis "y"

    set flag 0
    foreach keyitem [array names data] {
       if { [lsearch $keyitem tag ] != -1 } {
          set tag $data($keyitem)
          if { $data([ list $tag displayed ]) == 1 } {
               if { $data([ list $tag drawbars ])==1 } {
                   barplotter -tag $tag   ; # This plots a bar chart
               } elseif { $data([ list $tag is3d ]) == 0 } {
                   2dplotter -tag $tag    ; # This plots data points 
               } elseif { $data([ list $tag is3d ]) == 1 && $flag == 0 } {
                   3dplotter -tag $tag    ; # This plots free energy surfaces  
                   set flag 1
               } else { 
                   tk_messageBox -icon error -type ok -title Message -message "Can't plot more than one fes at once"   
                   return
               }
          }
       } 
    }
}

proc gtPlot::barplotter {args} {
    variable bounds
    variable data
    
    set tag [getargs $args "-tag" {}]

    # Get the number of points
    set np $data([ list $tag np ])

    # The base of the bars is the x-axis
    set base $bounds(yorigin_points)

    # Now draw the bar chart
    for { set i 0 } { $i<$np } { incr i } {
        set x1 [lindex $data([list $tag x]) $i]
        if { [expr $i + 1]==$np } { 
           set x2 [ expr 2*[lindex $data([list $tag x]) $i] - [lindex $data([list $tag x]) [expr $i - 1] ] ] 
        } else {
           set x2 [lindex $data([list $tag x]) [expr $i + 1] ]
        }
        set y [lindex $data([list $tag y]) $i]       
        drawLine -from $x1 $base -to $x1 $y -tag "$tag displayed bar" -fill $data([list $tag color]) -width 0.1 
        drawLine -from $x1 $y -to $x2 $y -tag "$tag displayed bar" -fill $data([list $tag color]) -width 0.1
        drawLine -from $x2 $y -to $x2 $base -tag "$tag displayed bar" -fill $data([list $tag color]) -width 0.1
   }
}

proc gtPlot::2dplotter {args} {
    variable bounds
    variable data
 
    set tag [getargs $args "-tag" {}]    

    # Get the number of points
    set np $data([ list $tag np ])

    # Draw data points
    for { set i 0 } { $i<$np } { incr i } {
        if { [expr $i%$data([list $tag stride])] == 0 } {
           set x [lindex $data([list $tag x]) $i]
           set y [lindex $data([list $tag y]) $i]
           if { $data([list $tag drawpoints ]) == 1 } {
              if { $x >= $bounds(xmin) && $x <= $bounds(xmax) && $y >= $bounds(ymin) && $y <= $bounds(ymax) } {
                   drawPixel -at $x $y -tag "pixel notfes $tag$i $tag displayed" -fill $data([list $tag color]) -outline black
              } else {
                   drawPixel -at $x $y -tag "pixel notfes $tag$i $tag" -fill $data([list $tag color]) -outline black
              }
           }
           if { $data([list $tag drawline]) == 1 && $i>0 } {
                set x0 [lindex $data([list $tag x]) [expr $i - 1] ]
                set y0 [lindex $data([list $tag y]) [expr $i - 1] ]
                drawLine -from $x0 $y0 -to $x $y -tag "$tag displayed dataline" -fill $data([list $tag color]) -width 0.1
           }
        }
    }

    if { $data([list $tag arehighlights]) == 1 } {
       drawHighlights -datatag $tag    ; # Extra work to deal with highlights
    }
}

proc gtPlot::3dplotter {args} {
    variable appearance
    variable bounds
    variable data
    variable canv
    variable ntics_d

    set tag [getargs $args "-tag" {}]

    # Check number of colors is specified properly
    if { ![ string is integer -strict $data([list $tag ncolors]) ] } { 
       tk_messageBox -icon error -type ok -title Message -message "Number of colors in scale not specified"    
       return
    }

    # Find the portion of the grid is inside the x and y range
    set gridbounds [ range2grid -tag $tag -axis x -range $bounds(xmin) $bounds(xmax) ]
    set xgridstart [lindex $gridbounds 0]   ; set xgridend [lindex $gridbounds 1]
    set gridbounds [ range2grid -tag $tag -axis y -range $bounds(ymin) $bounds(ymax) ]
    set ygridstart [lindex $gridbounds 0]   ; set ygridend [lindex $gridbounds 1]

    if { $xgridstart=="error" || $ygridstart=="error" } { return }

    # Find the range of z values
    if { $data([list $tag fes fixedz]) == 1 } {
       set bounds(zmin) $data([list $tag z min]) ; set bounds(zmax) $data([list $tag z max]) 
    } else { 
       set flag 0
       for { set i $xgridstart } { $i<$xgridend } { incr i } {
           for { set j $ygridstart } { $j<$ygridend } { incr j } {
               set dat1 [ lindex $data([list $tag fes $i $j]) 2 ]
               set dat2 [ lindex $data([list $tag fes $i [expr $j + 1] ]) 2 ]
               set dat3 [ lindex $data([list $tag fes [expr $i + 1] $j])  2 ]
               set dat4 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 2 ]
               set av [ expr ( $dat1 + $dat2 + $dat3 + $dat4 ) / 4. ]

               if { $flag == 0 } {
                   set bounds(zmin) $av ; set bounds(zmax) $av ; incr flag
               } elseif { $av > $bounds(zmax)} {
                   set bounds(zmax) $av
               } elseif { $av < $bounds(zmin)} { 
                   set bounds(zmin) $av
               }
           }
       }
       set appearance(zticspacing) [ format %1.0e [ expr ( $bounds(zmax) - $bounds(zmin) ) / $ntics_d ] ]
    }

    # Setup the colorscale
    set ncolors $data([list $tag ncolors])
    set colormap [ createColormap -colorscale $data([list $tag fescolors]) -ncolors $ncolors ]
    set color_incr [ expr ($bounds(zmax) - $bounds(zmin) ) / $ncolors ]
    createScalez -colorscale $data([list $tag fescolors]) -ncolors $ncolors  

    # Now plot the free energy surface
    for { set i $xgridstart } { $i<$xgridend } { incr i } {
        for { set j $ygridstart } { $j<$ygridend } { incr j } {
            # Find the point on the grid
            set x1 [ lindex $data([list $tag fes $i $j]) 0 ]
            set y1 [ lindex $data([list $tag fes $i $j]) 1 ]
            set x2 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 0 ]
            set y2 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 1 ]

            # Average the value of the fes at the four corners of this pixel
            set dat1 [ lindex $data([list $tag fes $i $j]) 2 ]
            set dat2 [ lindex $data([list $tag fes $i [expr $j + 1] ]) 2 ]
            set dat3 [ lindex $data([list $tag fes [expr $i + 1] $j])  2 ]
            set dat4 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 2 ]
            set av [ expr ( $dat1 + $dat2 + $dat3 + $dat4 ) / 4. ]

            # Convert this average value into a point on our discretized color scale
            set kk [ expr int( ( $av - $bounds(zmin) ) / $color_incr ) + 1 ]
            set color [ lindex $colormap [expr $kk - 1] ]

            # for { set kk 1 } { $kk<=$ncolors } { incr kk } {
            #     if { $av <= [ expr $zmin + $color_incr*$kk ] } { set color [ lindex $colormap [expr $kk - 1] ] ; break }
            # }
            # Draw this point
            drawFesPixel -topleft $x1 $y1 -botright $x2 $y2 -color $color -tag "pixel fes $tag"
        }
    }

    # Put the fes at the bottom of the stack
    $canv lower fes notfes 
}


proc gtPlot::range2grid {args} {
    variable data
    variable maxgrid

    set tag [getargs $args "-tag" {}]
    set axis [getargs $args "-axis" {}]
    set min [getargs $args "-range" {} 1]
    set max [getargs $args "-range" {} 2]

    if { $axis == "x" } {
        set n [expr $data([list $tag x ngrid]) - 1 ]   ; # Last point is skiped because of averages
        if { $n > $maxgrid } {
          tk_messageBox -icon error -type ok -title Message -message "Grid is too dense to plot in a reasonable time"
          return [list error]
        }
        set gmax $n ; set minf 0 
        for { set i 1 } { $i<$n } { incr i } {
           if { $minf==0 && [lindex $data([list $tag fes $i 0]) 0 ] >= $min } { set gmin [expr $i-1]; set minf 1 }
           if { $minf==1 && [lindex $data([list $tag fes $i 0]) 0] >= $max } { set gmax [expr $i]; break; }
        }
    } else {
        set n [expr $data([list $tag y ngrid]) - 1 ]  ; # Last point is skiped because of averages
        if { $n > $maxgrid } {
          tk_messageBox -icon error -type ok -title Message -message "Grid is too dense to plot in a reasonable time"
          return [list error]
        }
        set gmax $n ; set minf 0
        for { set i 1 } { $i<$n } { incr i } {
           if { $minf==0 && [lindex $data([list $tag fes 0 $i]) 1 ] >= $min } { set gmin [expr $i-1]; set minf 1 }
           if { $minf==1 && [lindex $data([list $tag fes 0 $i]) 1] >= $max } { set gmax [expr $i]; break; }
        }
    }

    return [list $gmin $gmax]
}

proc gtPlot::drawLine {args} {
   variable canv

   set x1 [getargs $args -from {} 1]
   set y1 [getargs $args -from {} 2]
   set x2 [getargs $args -to {} 1]
   set y2 [getargs $args -to {} 2]
   set color [getargs $args -fill {}]
   set tag [getargs $args "-tag" {}]
   set width [getargs $args "-width" {}]

   set x1 [ point2canvas -point $x1 -axis "x"]
   set y1 [ point2canvas -point $y1 -axis "y"]
   set x2 [ point2canvas -point $x2 -axis "x"]
   set y2 [ point2canvas -point $y2 -axis "y"]

   $canv create line $x1 $y1 $x2 $y2 -tag $tag -fill $color -width $width
}

proc gtPlot::drawFesPixel {args} {
   variable canv
   set x1 [getargs $args "-topleft" {} 1]  
   set y1 [getargs $args "-topleft" {} 2]
   set x2 [getargs $args "-botright" {} 1]
   set y2 [getargs $args "-botright" {} 2]
   set tag [getargs $args "-tag" {}]
   set color [getargs $args "-color" {}]
 
   # Tranform the input coordinates onto the canvas   
   set x1 [ point2canvas -point $x1 -axis "x"]
   set y1 [ point2canvas -point $y1 -axis "y"]
   set x2 [ point2canvas -point $x2 -axis "x"]
   set y2 [ point2canvas -point $y2 -axis "y"]

   $canv create rectangle $x1 $y1 $x2 $y2 -tag $tag -outline $color -fill $color
}


proc gtPlot::drawPixel {args} {
   variable canv
   variable appearance 
    
   set x [getargs $args "-at" {} 1]
   set y [getargs $args "-at" {} 2]
   set tag [getargs $args "-tag" {}]
   set fcolor [getargs $args "-fill" {}]
   set ocolor [getargs $args "-outline" {}]
 
   set x [ point2canvas -point $x -axis "x"]
   set y [ point2canvas -point $y -axis "y"]
   set s $appearance(size)
   
   $canv create rectangle [expr $x - $s] [expr $y - $s] [expr $x + $s] [expr $y + $s] -tag $tag -fill $fcolor -outline $ocolor
}

proc gtPlot::point2canvas {args} {
   variable appearance
   variable bounds
   variable canv
   
   set x [getargs $args "-point" {}]
   set d [getargs $args "-axis" {}]   
   
   if { $d == "x" } {
      set s [ expr [ winfo width $canv ] - 2.0*$appearance(pad) ]
      return [ expr $s * ( 0.5 + ( $x - $bounds(xcenter) ) / $bounds(xrange) ) + $appearance(pad) ]
   } else {
      set s [ expr [ winfo height $canv ] - 2.0*$appearance(pad) ] 
      return [ expr $s * ( 0.5 - ( $x - $bounds(ycenter) ) / $bounds(yrange) ) + $appearance(pad) ]     
   }   
}

proc gtPlot::canvas2point {args} {
   variable appearance
   variable bounds
   variable canv
   
   set x [getargs $args "-point" {}]
   set d [getargs $args "-axis" {}]   
   
   if { $d == "x" } {
      set s [ expr [ winfo width $canv ] - 2.0*$appearance(pad) ]
      return [ expr ( ( $x - $appearance(pad) ) / $s - 0.5 ) * $bounds(xrange) + $bounds(xcenter) ]
   } else {
      set s [ expr [ winfo height $canv ] - 2.0*$appearance(pad) ] 
      return [ expr ( 1.0 - ( $x - $appearance(pad) ) / $s - 0.5 ) * $bounds(yrange) + $bounds(ycenter) ] 
   }   
}

proc gtPlot::FesStart {canv x y} {
   variable fesID

   set st $fesID(st)
   if { $st == 0 } {
     $canv delete fesint
     set fesID([expr $st+0]) [$canv create rect $x $y $x $y -tag fesint]
   } else {
     set fesID([expr $st+0]) [$canv create rect $x $y $x $y -dash { 4 4 } -tag fesint]
   }
   set fesID([expr $st+1]) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag "feshand fesint"]
   set fesID([expr $st+2]) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag "feshand fesint"]
   set fesID(xmin) $x     ;   set fesID(xmax) $x
   set fesID(ymin) $y     ;   set fesID(ymax) $y
}

proc gtPlot::FesMove {canv x y} {
   variable fesID 
   
   set st $fesID(st)
   $canv coords $fesID([expr $st+0]) [lreplace [$canv coords $fesID([expr $st+0])] 2 3 $x $y]
   $canv coords $fesID([expr $st+2]) [lreplace [$canv coords $fesID([expr $st+2])] 0 3 [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] ]
   set fesID(xmax) $x    ;    set fesID(ymin) $y
}

proc gtPlot::FesEnd {canv x y} {
   variable fesID
   variable data

   $canv delete feshand

   # Keep track of whether this is the first or second box drawn
   if { $fesID(st)==0 } {
     set fesID(st) 3
     # Save the first set of xmax and xmin data
     set fesID(sxmin) $fesID(xmin) ; set fesID(sxmax) $fesID(xmax)
     set fesID(symin) $fesID(ymin) ; set fesID(symax) $fesID(ymax)
   } else {
     set fesID(st) 0
   }

   writeFesData $canv     ; # Actually do the integration of the fes
}

proc gtPlot::writeFesData {canv} {
   variable appearance
   variable fesID
   variable data
   variable kT

   # Check that kT is set
   if { $kT(kT)==0 } {
      $canv delete fesint
      ::gtPlot::setkT
      return
   }

   # Check that there are some rectangles to integrate inside
   if { [llength [$canv find withtag fesint] ]==0 } { return }

   # Find the fes we have to integrate
   set flag 0
   foreach keyitem [array names data] {
      if { [lsearch $keyitem tag ] != -1 } {
         set tag $data($keyitem)
         if { $data([ list $tag displayed ]) == 1 } {
              if { $data([ list $tag isfes ]) == 1 } {
                 set flag $tag
              } elseif { $data([ list $tag isfes ]) == 1 && $flag!=0 } {
                 $canv delete fesint
                 tk_messageBox -icon error -type ok -title Message -message "More than one fes plotted - I don't know what to do"
                 return
              }
         }
      }
   }
   if { $flag==0 } { $canv delete fesint ; return }

   # Delete the old fes data from the canvas
   $canv delete festext 

   # Outer loop here decides how many integrations we have to do
   if { $fesID(st)==3 } { 
     # Do the integration of fes
     if { $data([list $tag is3d])==1 } {
         set fes [ integrate3dFes -datatag $flag -xrange $fesID(sxmin) $fesID(sxmax) -yrange $fesID(symin) $fesID(symax) ]
     } else {
         set fes [ integrate2dFes -datatag $flag -xrange $fesID(sxmin) $fesID(sxmax) ]
     }
     # Write on the canvas the value of the fes
     $canv create text 10 [ expr [ winfo height $canv ] - 20 ] -font $appearance(font) -anchor sw -text "G(solid box) = [format %.4g $fes] \n\n" -tag "festext fesint"
   } else {
     # Do the integration of the two regions of the fes that need integrating 
     if { $data([list $tag is3d])==1 } {
        set fes1 [ integrate3dFes -datatag $flag -xrange $fesID(sxmin) $fesID(sxmax) -yrange $fesID(symin) $fesID(symax) ]
        set fes [ integrate3dFes -datatag $flag -xrange $fesID(xmin) $fesID(xmax) -yrange $fesID(ymin) $fesID(ymax) ]
     } else {
        set fes1 [ integrate2dFes -datatag $flag -xrange $fesID(sxmin) $fesID(sxmax) ]
        set fes [ integrate2dFes -datatag $flag -xrange $fesID(xmin) $fesID(xmax) ]
     }
     # Write on the canvas
     $canv create text 10 [ expr [ winfo height $canv ] - 20 ] -font $appearance(font) -anchor sw -text " G(solid box) = [format %.4g $fes1] \n G(dashed box) = [format %.4g $fes] \n difference = [format %.4g [expr $fes1-$fes]]" -tag "festext fesint" -justify left
   }
}

proc gtPlot::integrate2dFes {args} {
   variable data
   variable kT

   # We must get the datatag 
   set tag [getargs $args "-datatag" {}]
   set rxmin [getargs $args "-xrange" {none} 1]
   set rxmax [getargs $args "-xrange" {none} 2]

   # This gets the coordinates of the fes box
   set xmin [canvas2point -point $rxmin -axis x ]
   set xmax [canvas2point -point $rxmax -axis x ]

   # Now integrate the fes in the highlighted region
   set v 0 ; set np $data([ list $tag np ])
   for {set i 0} { $i<$np } { incr i } {
      set x [lindex $data([list $tag x]) $i]
      if { $x>$xmin && $x<$xmax } {
         set z [lindex $data([list $tag y]) $i]
         set v [ expr $v + exp(-$z/$kT(kT)) ]
      }
   } 
   return [ expr -$kT(kT)*log($v) ]
}


proc gtPlot::integrate3dFes {args} {
   variable data
   variable kT

   # We must get the datatag 
   set tag [getargs $args "-datatag" {}]
   set rxmin [getargs $args "-xrange" {none} 1]
   set rxmax [getargs $args "-xrange" {none} 2]
   set rymin [getargs $args "-yrange" {none} 1]
   set rymax [getargs $args "-yrange" {none} 2]

   # This gets the coordinates of the fes box
   set xmin [canvas2point -point $rxmin -axis x ]
   set ymin [canvas2point -point $rymin -axis y ]
   set xmax [canvas2point -point $rxmax -axis x ]
   set ymax [canvas2point -point $rymax -axis y ]

   # Convert the size of the box to something on the grid
   set gridbounds [ range2grid -tag $tag -axis x -range $xmin $xmax ]
   set xgridstart [lindex $gridbounds 0]   ; set xgridend [lindex $gridbounds 1]
   set gridbounds [ range2grid -tag $tag -axis y -range $ymin $ymax ]
   set ygridstart [lindex $gridbounds 0]   ; set ygridend [lindex $gridbounds 1]

   if { $xgridstart=="error" || $ygridstart=="error" } { return }

   # Now integrate the fes in the highlighted region
   set v 0
   for { set i $xgridstart } { $i<$xgridend } { incr i } {
        for { set j $ygridstart } { $j<$ygridend } { incr j } {
           set dat1 [ lindex $data([list $tag fes $i $j]) 2 ]
           set dat2 [ lindex $data([list $tag fes $i [expr $j + 1] ]) 2 ]
           set dat3 [ lindex $data([list $tag fes [expr $i + 1] $j])  2 ]
           set dat4 [ lindex $data([list $tag fes [expr $i + 1] [expr $j + 1] ]) 2 ]
           set av [ expr ( $dat1 + $dat2 + $dat3 + $dat4 ) / 4. ]
           set v [ expr $v + exp(-$av/$kT(kT)) ]
        }
   }
   return [ expr -$kT(kT)*log($v) ]
}

proc ::gtPlot::calcFes {args} {
   variable fesID
   variable data
   variable canv
   variable kT

   set festag [getargs $args "-festag" {}]
   set dtag [getargs $args "-datatag" {}]

   # Check that kT is set
   if { $kT(kT)==0 } {
      $canv delete fesint
      ::gtPlot::setkT
      return "error"
   }

   # Check that the correct fes is plotted   
   if { $data([ list $festag displayed ]) != 1 } {
        tk_messageBox -icon error -type ok -title Message -message "The fes plotted does not match the one desired for the difference"
        return "error"
   }
   
   # Now check that there is a highlighted region
   if { [llength [$canv find withtag fesint] ]==0 } {
      tk_messageBox -icon error -type ok -title Message -message "You must have highlighted some regions of the free energy surface using cntr-lclick to examine convergence"
      return "error"
   }  

   # This will do the integration of the fes (in box 1
   if { $data([list $dtag is3d])==1 } {
      set fes1 [ integrate3dFes -datatag $dtag -xrange $fesID(sxmin) $fesID(sxmax) -yrange $fesID(symin) $fesID(symax) ]
   } else {
      set fes1 [ integrate2dFes -datatag $dtag -xrange $fesID(sxmin) $fesID(sxmax) ]
   }
   
   if { $fesID(st)==3 } { return $fes1 }

   # This will do the integration of the fes (in box 2
   if { $data([list $dtag is3d])==1 } {
      set fes2 [ integrate3dFes -datatag $dtag -xrange $fesID(xmin) $fesID(xmax) -yrange $fesID(ymin) $fesID(ymax) ]
   } else {
      set fes2 [ integrate2dFes -datatag $dtag -xrange $fesID(xmin) $fesID(xmax) ]
   }

   return [expr $fes1-$fes2]
}

proc gtPlot::FesClear {canv x y} {
   variable fesID
   $canv delete fesint 
   set fesID(st) 0
}

proc gtPlot::ZoomStart {canv x y} {
   variable zoomID
   set zoomID(0) [$canv create rect $x $y $x $y -dash { 4 4 } -tag zoom ]
   set zoomID(1) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag zoom ]
   set zoomID(2) [$canv create rect [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] -tag zoom ]
   set zoomID(xmin) $x     ;   set zoomID(xmax) $x
   set zoomID(ymin) $y     ;   set zoomID(ymax) $y
}

proc gtPlot::ZoomMove {canv x y} {
   variable zoomID
   $canv coords $zoomID(0) [lreplace [$canv coords $zoomID(0)] 2 3 $x $y]
   $canv coords $zoomID(2) [lreplace [$canv coords $zoomID(2)] 0 3 [expr $x-2] [expr $y-2] [expr $x+2] [expr $y+2] ] 
   set zoomID(xmax) $x    ;    set zoomID(ymin) $y 
}

proc gtPlot::ZoomEnd {canv x y} {
   variable zoomID
   variable data

   # Delete the zoom box
   $canv delete zoom
   $canv delete pixel 
   $canv delete axis
   $canv delete whereIam
   $canv delete highlight_line   
   $canv delete basinsize  
   $canv delete dataline  

   # A check that there is some data plotted
   set flag 0
   foreach keyitem [array names data] {
      # Look for all the tag names
      if { [ lsearch $keyitem tag ] != -1 } {
           set tag $data($keyitem) 
           # Check if data is displayed
           if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
      }
   }
   if { $flag == 0 } { return } 

   # This gets the coordinates of the zoom box
   set xmin [canvas2point -point $zoomID(xmin) -axis x ]
   set ymin [canvas2point -point $zoomID(ymin) -axis y ]  
   set xmax [canvas2point -point $zoomID(xmax) -axis x ]   
   set ymax [canvas2point -point $zoomID(ymax) -axis y ]

   # Reset the axis
   setAxis -xrange $xmin $xmax -yrange $ymin $ymax

   # Redraw the axis
   drawAxis -axis "x"
   drawAxis -axis "y"

   # This plots all the data sets
   set flag 0
   foreach keyitem [array names data] {
      if { [lsearch $keyitem tag ] != -1 } {
         set tag $data($keyitem)
         if { $data([ list $tag displayed ]) == 1 } {
            if { $data([ list $tag drawbars ])==1 } {
                barplotter -tag $tag   ; # This plots a bar chart
            } elseif { $data([ list $tag is3d ]) == 0 } {
                2dplotter -tag $tag    ; # This plots data points 
            } elseif { $data([ list $tag is3d ]) == 1 && $flag == 0 } {
                3dplotter -tag $tag    ; # This plots free energy surfaces  
                set flag 1
            } else {
                tk_messageBox -icon error -type ok -title Message -message "Can't plot more than one fes at once"   
                return
            }
         }
      }
   }
}

proc gtPlot::writeCoord {canv x y} {
   variable appearance   
   variable data

   # Delete the old label   
   $canv delete whereIam

   # A check that there is some data plotted
   set flag 0
   foreach keyitem [array names data] {
      # Look for all the tag names
      if { [ lsearch $keyitem tag ] != -1 } {
           set tag $data($keyitem)
           # Check if data is displayed
           if { $data([ list $tag displayed ]) == 1 } { set flag 1 }
      }
   }
   if { $flag == 0 } { return }

   # This gets our coordinates
   set px [format %.4g [canvas2point -point $x -axis x ] ]
   set py [format %.4g [canvas2point -point $y -axis y ] ]

   # Write the position on the canvas
   $canv create text [ expr [ winfo width $canv ] - 5] 15 -font $appearance(font) -anchor se -text "($px, $py)" -tag "notfes whereIam"  
}

proc ::gtPlot::printCanvas {args} {
  variable canv
        
  set filename [tk_getSaveFile \
        -initialfile "colvars.ps" \
        -title "Print graph to file" \
        -parent $canv \
        -filetypes [list {{Postscript files} {.ps}} {{All files} {*}}]]
  if {$filename != ""} {
    $canv postscript -file $filename
  }
  return

}

# This sets up our colormaps
proc gtPlot::createColormap {args} {
  set colorscale [getargs $args "-colorscale" {}]
  set ncolors [getargs $args "-ncolors" {}] 

  switch -- $colorscale {
    hsv {
        set hueStart     0.0
        set hueEnd     240.0
        set colorMap   {}
    
        for {set i 0} {$i <= $ncolors} {incr i} {
            set dh [expr {($hueStart - $hueEnd) / ($ncolors - 1)}]
            set hue  [expr {$hueStart - ($i * $dh)}]
            if {$hue < 0.0} {
                set hue  [expr {360.0 + $hue}]
            }
            set rgbList [Hsv2rgb $hue 1.0 1.0]
            set r    [expr {int([lindex $rgbList 0] * 65535)}]
            set g    [expr {int([lindex $rgbList 1] * 65535)}]
            set b    [expr {int([lindex $rgbList 2] * 65535)}]
    
            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    hot {
        set colorMap {}
        set nc1          [expr {int($ncolors * 0.33)}]
        set nc2          [expr {int($ncolors * 0.67)}]
    
        for {set i 0} {$i <= $ncolors} {incr i} {
    
            if {$i <= $nc1} {
                set fval  [expr { double($i) / (double($nc1)) } ]
                set r     [expr {int($fval * 65535)}]
                set g     0
                set b     0
            } else {
                if {$i <= $nc2} {
                    set fval  [expr { double($i-$nc1) / (double($nc2-$nc1)) } ]
                    set r     65535
                    set g     [expr {int($fval * 65535)}]
                    set b     0
                } else {
                    set fval  [expr { double($i-$nc2) / (double($ncolors-$nc2)) } ]
                    set r     65535
                    set g     65535
                    set b     [expr {int($fval * 65535)}]
                }
            }
            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    cool {
        set colorMap {}
        for {set i 0} {$i <= $ncolors} {incr i} {
    
            set fval  [expr { double($i) / (double($ncolors)-1) } ]
            set val   [expr { 1.0 - $fval }]
    
            set r    [expr {int($fval * 65535)}]
            set g    [expr {int($val * 65535)}]
            set b    65535
    
            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    grey {  
        set colorMap {}
        for {set i 0} {$i <= $ncolors} {incr i} {

            set fval  [expr { double($i) / (double($ncolors)-1) } ]
            set val  [expr {0.4 + (0.5 * $fval) }]

            set r    [expr {int($val * 65535)}]
            set g    [expr {int($val * 65535)}]
            set b    [expr {int($val * 65535)}]

            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
        }
        return $colorMap
    }
    jet {
         set hueStart   240.0
         set hueEnd       0.0
         set colorMap   {}

         for {set i 0} {$i <= $ncolors} {incr i} {

            set dh [expr {($hueStart - $hueEnd) / ($ncolors - 1)}]
            set hue  [expr {$hueStart - ($i * $dh)}]
            if {$hue < 0.0} {
                set hue  [expr {360.0 + $hue}]
            }
            set rgbList [Hsv2rgb $hue 1.0 1.0]
            set r    [expr {int([lindex $rgbList 0] * 65535)}]
            set g    [expr {int([lindex $rgbList 1] * 65535)}]
            set b    [expr {int([lindex $rgbList 2] * 65535)}]

            set color  [format "#%.4x%.4x%.4x" $r $g $b]
            lappend colorMap $color
         } 
         return $colorMap 
    }
  }
}

proc gtPlot::Hsv2rgb {h s v} {
    set v [expr {double($v)}]
    set r [set g [set b 0.0]]
    if {$h == 360} { set h 0 }
    # if you feed the output of rgb2hsv back into this
    # converter, h could have the value -1 for
    # grayscale colors.  Set it to any value in the
    # valid range.
    if {$h == -1} { set h 0 }
    set h [expr {$h/60}]
    set i [expr {int(floor($h))}]
    set f [expr {$h - $i}]
    set p1 [expr {$v*(1-$s)}]
    set p2 [expr {$v*(1-($s*$f))}]
    set p3 [expr {$v*(1-($s*(1-$f)))}]
    switch -- $i {
        0 { set r $v  ; set g $p3 ; set b $p1 }
        1 { set r $p2 ; set g $v  ; set b $p1 }
        2 { set r $p1 ; set g $v  ; set b $p3 }
        3 { set r $p1 ; set g $p2 ; set b $v  }
        4 { set r $p3 ; set g $p1 ; set b $v  }
        5 { set r $v  ; set g $p1 ; set b $p2 }
    }
    return [list $r $g $b]
}
