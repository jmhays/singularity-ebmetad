#
# Graphical user interface for PLUMED
#
# $Id: plumed_gui.tcl,v 1.0 2010/11/15
# 

# To do list
# (1) Need to do a lot of testing of using gui in tandem with multiple molecules in vmd (this needs to be more robust methinks)
# (2) Must add keywords for all known colvars to cv interpreter    --    PCA

package provide plumed_gui 1.0
package require gtPlot 1.0

vmd_install_extension plumed_gui plumedVis_tk "Analysis/Plumed GUI"

proc plumedVis_tk {} {
    plumedVis::create_gui
    return $plumedVis::w
}

namespace eval plumedVis:: {
  namespace export CVPlot 
  variable molid      -1        ;# We store the current molecule that is on top here
  variable w                    ;# handle to window
  variable cv                   ;# handle to cv appearance window
  variable fes                  ;# handle to fes appearance window
  variable filename             ;# filename to read in
  variable cvstride     1       ;# stride between displayed cv values
  variable highlight_nbefore 0  ;# Number of frames to highlight before the current point
  variable highlight_nafter  0  ;# Number of frames to highlight after the current point 
  variable highlight_lines   1  ;# Draw lines between data points in highlight
  variable xcoord               ;# The coord on the x axis
  variable ycoord               ;# The coord on the y axis
  variable cvdata               ;# an array containing all the CV data in lists
  variable recondata            ;# an array containing all the basin data in lists
  variable fesdata              ;# an array containing all fesdata plotted
  variable appearance           ;# An array of stuff that controls the appearance
  variable readcolvars          ;# A binary flag that tells us if we have colvar info
  variable allcolvars           ;# A binary flag that tells us if we are plotting all colvars or just one
  variable sumhills             ;# This is an array that holds the data we need to sum hills
  variable lenunit              ;# The unit of length we are using
  variable tfilename            ;# File containing masses and charges
  variable recon_init_size  3   ;# This is the initial size of reconnaissance basins (i.e. SIZE_PARAM)
}

proc plumedVis::setDefaults {args} {
   variable appearance
   variable fesdata
   variable readcolvars
   set appearance(canvasSize)     400           ; #The size of the canvas
   set appearance(canvasPadding)  60            ; #The ammount of white space around the canvas
   set appearance(boxScale)       0.25          ; #This is the scale for box sizes
   set appearance(font)          "Helvetica 8"  ; #Font for axis labels and so on 
   # These are things that can be adjusted
   set appearance(trajcol)       "grey"         ; #This is the colorscale used for drawing the trajectory (grey/jet)
   set appearance(nfescolors)     100           ; #This is the number of colors used to draw free energy surfaces
   set appearance(fescol)        "jet"          ; #This is the colorscale used for drawing the fes (grey/jet)
   set appearance(boxSize)        10            ; #Default size of boxes
   set appearance(fixz)            0            ; #Dont change z axis when we zoom and so on
   set appearance(plotdiff)        0            ; #Plot the difference between free energy surfaces

   set readcolvars      0   ; # Currently no colvars have been read in   (need to make this work with multiple molecules in vmd)        GAT
   set allcolvars       0   ; # This keeps track of whether we are in all cv plotting mode or just in regular (2cv) plotting mode
}

proc plumedVis::create_gui {} {
  # Just create the window and initialize data structures
  # No molecule has been selected yet
  # Also set up traces on VMD variables

  variable appearance
  variable w
  variable molid
  variable filename
  variable xcoord
  variable ycoord
  variable fesdata

  global vmd_frame
  global env

  # Check plumedir is defined (if this is not the case we should crash
  if { $env(plumedir) == "" } { tk_messageBox -icon error -type ok -title Message -message "I have not found plumedir in env update .bashrc" }

  # If already initialized, just turn on 
  if [winfo exists .bplot] {
    wm deiconify $w
    return
  }

  set w [toplevel ".bplot"]
  wm title $w "CVPlot - Dynamic collective variable plots for VMD"
  wm resizable $w 1 1
  bind $w <Destroy> plumedVis::destroy

  # Set all the defaults
  setDefaults

  # Create a tempory directory in which to do calculations 
  set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"
  file mkdir $tmpd 

  frame $w.top -padx 1m -pady 1m

  pack [ createMenubar $w.top.menubar ] -padx 1 -fill x -side top

  # This is our graph object
  ::gtPlot::create $w.fr -size $appearance(canvasSize) -font $appearance(font) -padding $appearance(canvasPadding) 
  pack $w.fr -in $w.top -fill both -expand true -side top
  # Add axis controls to appearance menu   
  $w.top.menubar.look.menu add checkbutton -label "Fix z scale" -variable "plumedVis::appearance(fixz)"
  $w.top.menubar.look.menu add checkbutton -label "Plot fes difference" -variable "plumedVis::appearance(plotdiff)" 
  $w.top.menubar.look.menu add command -label "Automatic axis" -command [namespace code ::gtPlot::autoAxes ]
  $w.top.menubar.look.menu add command -label "Axis properties" -command [namespace code ::gtPlot::axisPropertiesDialogue ]

  # This intializes the spinbox menus
  lappend fesdata(frames) 1 ;   set fesdata(fesframe) 1 

  frame $w.top.pmenus -relief raised -bd 2
  pack [ axisMenu $w.top.pmenus.xaxis -text "x-axis" -textvariable plumedVis::xcoord ] -in $w.top.pmenus -side left
  pack [ axisMenu $w.top.pmenus.yaxis -text "y-axis" -textvariable plumedVis::ycoord ] -in $w.top.pmenus -side left
  pack [ button $w.top.pmenus.allcvs -text "all cvs" -relief raised -command plumedVis::plotAllCVs ] -in $w.top.pmenus -side left  
  pack [ label $w.top.pmenus.slab -text "fes frame" -anchor e ] -in $w.top.pmenus -side left
  eval spinbox $w.top.pmenus.spinb -textvariable plumedVis::fesdata(fesframe) -values $plumedVis::fesdata(frames) -width 5
  pack $w.top.pmenus.spinb -in $w.top.pmenus -side left 
  pack $w.top.pmenus -in $w.top -padx 1 -fill x -side top
  pack $w.top -fill both -expand true                      ; # Complete the packing

  # This binds every box to a frame
  $w.fr.canvas bind colvar: <Button-1> [namespace code { plumedVis::Goto %W } ]
  # This binds every basin to the size drawer (in theory)
  $w.fr.canvas bind basins: <Button-1> [namespace code { plumedVis::DrawBasinSize %W } ]

  # This controls the plotter
  trace add variable plumedVis::xcoord write [namespace code plotData]
  trace add variable plumedVis::ycoord write [namespace code plotData]
  trace add variable plumedVis::appearance(fixz) write [namespace code setupZoom]
  # Update the highlight which indicates where we are in the trajectory
  trace add variable vmd_frame write [namespace code updateHighlight]

  set molid [ findMolecule ]    ; # Find the molecule on top 
}

proc plumedVis::DrawBasinSize {w} {
  variable molid
  variable recondata

  set id [ $w find withtag current]
  set taglist [$w gettags $id]
  set listindex [lsearch -glob $taglist basins*]
  if { $listindex < 0 } { return }
  set frametag [lindex $taglist $listindex]
  lassign [split $frametag :] foo bas

  drawBasin -basin $bas
  if { $recondata([list $molid basin_molid]) == $molid } { animate goto $bas }
}

proc plumedVis::drawBasin {args} {
  variable recondata
  variable cvdata
  variable molid
  variable xcoord
  variable ycoord

  set bas [getargs $args "-basin" {}]
  set basmol $recondata([list $molid basin_molid])

  set center [ list [lindex $cvdata([list $basmol $xcoord]) $bas] [lindex $cvdata([list $basmol $ycoord]) $bas] ]
  set shape [list [lindex $cvdata([list shape $basmol $xcoord $xcoord]) $bas] [lindex $cvdata([list shape $basmol $xcoord $ycoord]) $bas] \
                  [lindex $cvdata([list shape $basmol $ycoord $xcoord]) $bas] [lindex $cvdata([list shape $basmol $ycoord $ycoord]) $bas] ]

  ::gtPlot::drawBasinSize -center $center -shape $shape -size $cvdata([list $basmol $bas sizes]) 
}

proc plumedVis::destroy {args} {
   variable xcoord
   variable ycoord

   global vmd_frame

   # Remove the traces
   trace remove variable plumedVis::xcoord write [namespace code plotData]
   trace remove variable plumedVis::ycoord write [namespace code plotData]
   trace remove variable plumedVis::appearance(fixz) write [namespace code setupZoom]
   trace remove variable plumedVis::fesdata(fesframe) write [namespace code updateFes]
   trace remove variable vmd_frame write [namespace code updateHighlight]
}

proc plumedVis::Goto {w} {
  set id [ $w find withtag current]
  set taglist [$w gettags $id]
  set listindex [lsearch -glob $taglist colvar*]
  if { $listindex < 0 } { return }
  set frametag [lindex $taglist $listindex]
  lassign [split $frametag :] foo frame
  animate goto $frame
}

proc plumedVis::updateHighlight {args} {
  variable w
  variable molid
  variable cvdata
  variable recondata
  variable appearance
  variable highlight_nbefore
  variable highlight_nafter
  variable highlight_lines 
  variable allcolvars

  global vmd_frame

  # Create new menus if molecule has changed (also deletes the plot 
  if { $molid != [ findMolecule ] } { set molid [findMolecule] ; updateMenus ; return }

  # Check for non insane values of highlight_nbefore and highlight_nafter
  if { ![string is integer -strict $highlight_nbefore] || $highlight_nbefore<0 } { return }
  if { ![string is integer -strict $highlight_nafter] || $highlight_nafter<0 } { return }

  # If we are in allcv plotting mode then we just redraw the whole bar chart.
  if { $allcolvars==1 } { plotAllCVs ; return } 

  # This actually draws the highlight
  ::gtPlot::updateHighlight -datatag colvar: -highlight $vmd_frame($molid) -lines $highlight_lines -nbefore $highlight_nbefore -nafter $highlight_nafter -colors $appearance(trajcol)

  # Check for recondata
  set flag 0
  foreach cvname $cvdata([list $molid cv_names]) {
      set flag [expr $flag + $recondata([ list $molid $cvname isdata]) ] 
  }  

  # This makes sure the basins always appear on top -- GAT maybe put a logical around this 
  if { $flag != 0 } {
     $w.fr.canvas raise basins: colvar: 
     $w.fr.canvas raise highlighted basins:

     # Test if we are basins as we would like to draw sizes if we are 
     set basmol $recondata([list $molid basin_molid])
     if { $basmol==$molid } { drawBasin -basin $vmd_frame($molid) }
  } 
}

# This plots a bar chart of the cvs in the current frame
proc plumedVis::plotAllCVs {args} {
   variable fesdata
   variable readcolvars
   variable allcolvars
   variable xcoord
   variable ycoord
   variable cvdata
   variable molid

   global vmd_frame

   # Get rid of any plots of cvs in 2d
   if { $xcoord!="none" } { set xcoord "none" }
   if { $ycoord!="none" } { set ycoord "none" }
   # Get rid of any free energy surface from the plot
   ::gtPlot::remove -datatag free_energy

   # Now make sure there is some colvar data to plot
   if { $readcolvars==0 } { return }

   # Check that molid is on top
   if { $molid != [ findMolecule ] } { set molid [findMolecule] }

   # Extract the colvar data for this particular frame
   set n 0 
   foreach cvname [lsearch -all -inline $cvdata([ list $molid cv_names ]) cv* ] {
      lappend xdata $n 
      lappend ydata [lindex $cvdata([list $molid $cvname]) $vmd_frame($molid)]
      incr n
   } 

   # Remove the old allcolvars data
   ::gtPlot::remove -datatag allcolvars
   # Pass data to gtplot
   ::gtPlot::readData -datatag allcolvars -npoints [llength $ydata] -xdata $xdata -ydata $ydata 
   # Draw the graph
   ::gtPlot::plot2d -datatag allcolvars -stride 1 -color black -drawbars

   set allcolvars 1   ;  # This tells all other routines that all the colvars are plotted
}

# This plots the data 
proc plumedVis::plotData {args} {
   variable molid
   variable appearance
   variable cvdata
   variable recondata
   variable cvstride
   variable highlight_nbefore
   variable highlight_nafter
   variable highlight_lines
   variable xcoord
   variable ycoord
   variable allcolvars

   global vmd_frame

   # Delete the bar charts of all the colvars 
   ::gtPlot::remove -datatag allcolvars ; set allcolvars 0  

   # Delete data and return if nothing is set
   if { $xcoord=="none" || $ycoord=="none" } { 
        ::gtPlot::remove -datatag colvar: 
        ::gtPlot::remove -datatag basins:
        return
   }

   # Check that molid is on top
   if { $molid != [ findMolecule ] } { set molid [findMolecule] ; updateMenus ; return }

   # Check that data is available for this molecule
   set flag 0
   foreach keyitem [array names cvdata] {
      if { [lsearch $keyitem $molid ] != -1 && [lsearch $keyitem $xcoord] != -1 && [lsearch $keyitem shape] == -1 } { incr flag }
      if { [lsearch $keyitem $molid ] != -1 && [lsearch $keyitem $ycoord] != -1 && [lsearch $keyitem shape] == -1 } { incr flag } 
   } 

   # If no data available purge the menus and return  
   if { $flag!=2 } { updateMenus; tk_messageBox -icon error -type ok -title Message -message "No cv information available for this molecule"; return; }

   # Pass the data to the graph plotter
   ::gtPlot::readData -datatag colvar: -npoints [molinfo $molid get numframes] -xdata $cvdata([list $molid $xcoord]) -ydata $cvdata([list $molid $ycoord])

   # Set the highlight in gtPlot (this is needed for zoom and so on)
   ::gtPlot::setHighlights -datatag colvar: -highlight $vmd_frame($molid) -lines $highlight_lines -nbefore $highlight_nbefore -nafter $highlight_nafter -colors $appearance(trajcol)

   # Plot the data - highlighting is also done by this routine as we have set this up above
   ::gtPlot::plot2d -datatag colvar: -pointsize [expr $appearance(boxScale)*$appearance(boxSize)] -stride $cvstride -color white -drawpoints

   # If there is reconnaissance data for this molecule and coords plot it
   if { $recondata([list $molid $xcoord isdata])==1 && $recondata([list $molid $ycoord isdata])==1 } {
      set basmol $recondata([list $molid basin_molid])
      ::gtPlot::readData -datatag basins: -npoints [molinfo $basmol get numframes] -xdata $cvdata([list $basmol $xcoord]) -ydata $cvdata([list $basmol $ycoord])
      ::gtPlot::plot2d -datatag basins: -pointsize [expr $appearance(boxScale)*$appearance(boxSize)] -stride 1 -color blue -drawpoints
   }
}

proc plumedVis::storeSelection {args} {
   variable w
   variable molid
   variable recondata
   variable xcoord
   variable ycoord
   variable cvdata
  
   global vmd_frame

   # Return if nothing is set
   if { $xcoord=="none" || $ycoord=="none" } { return }

   # Find the tmp directory in which to do i/o
   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"
   set owd [pwd]
   cd $tmpd

   # Now find the molecules that are in the selected area
   set nfr 0
   foreach id [ $w.fr.canvas find withtag displayed ] {
     # This converts the tags of the molecules in the selection
     # To a frame number
     set taglist [$w.fr.canvas gettags $id]
     if { [lsearch $taglist colvar*]>=0  } {
        set listindex [lsearch -glob $taglist colvar*]
        if { $listindex < 0 } { return }
        set frametag [lindex $taglist $listindex]
        lassign [split $frametag :] foo frame

        # This will output coordinate information to a pdb file
        [ atomselect $molid all frame $frame ] writepdb select_tmp_$nfr.pdb 

        incr nfr  ;# Keep track of the number of frames in the selection
     }
   }

   if { $nfr > 0 } { 
      # Turn of the trace on vmd_frame to prevent plotData being called during readin of new files
      trace remove variable vmd_frame write [namespace code updateHighlight]

      # tell vmd to read in the data and load it as a new molecule (N.B. the selection will be the new top molecule)
      mol new select_tmp_0.pdb type pdb 
      for { set i 1 } { $i<$nfr } { incr i } { mol addfile select_tmp_$i.pdb type pdb }

      set newmol [findMolecule]   ; # Get the molid of the new molecule (for tagging data)
      mol delrep 0 $newmol        ; # Delete the default representation of the new molecule (we copy those of the old)

      #puts "Check in abstract 0 [ molinfo $molid get { {rep 0} {selection 0} {color 0} } ]"
      #puts "Check in abstract 1 [ molinfo $molid get { {rep 1} {selection 1} {color 1} } ]"

      # Get the list of representations for the old molecule 
      set nreps [molinfo $molid get numreps]
      for { set i 0 } { $i<$nreps } { incr i } {
          set toget "{ rep $i} {selection $i} {color $i}"
          set rep [ molinfo $molid get $toget ]
          mol addrep $newmol 
          mol modstyle $i $newmol [lindex $rep 0]
          mol modselect $i $newmol [lindex $rep 1]
          mol modcolor $i $newmol [lindex $rep 2]
      } 
 
      # Copy the cvdata of the old molecule to the new molecule (only relevant frames)
      set cvdata([list $newmol cv_names]) $cvdata([list $molid cv_names])
      foreach id [ $w.fr.canvas find withtag displayed ] {
         set taglist [$w.fr.canvas gettags $id]
         if { [lsearch $taglist colvar*]>=0  } {
            # This converts the tags of the molecules in the selection
            # To a frame number
            set listindex [lsearch -glob $taglist colvar*]
            if { $listindex < 0 } { return }
            set frametag [lindex $taglist $listindex]
            lassign [split $frametag :] foo frame

            foreach cv_name $cvdata([ list $molid cv_names ]) {
               lappend cvdata([ list $newmol $cv_name ]) [lindex $cvdata([list $molid $cv_name]) $frame] 
            }
         }
      }
      # There is no recondata for the new molecule
      foreach cvname $cvdata([list $newmol cv_names]) { set recondata([ list $newmol $cvname isdata]) 0 }

      set molid $newmol    ;# The top molecule has changed on read in
      plotData             ;# We replot the trajectory for the new molecule 

      # Turn the trace back on
      trace add variable vmd_frame write [namespace code updateHighlight]

      # Clean up all the pdb files we have created  
      eval file delete [glob -dir [pwd] select_tmp_*] 
   }

   cd $owd    ; # Go back to the original directory
}

# This finds the current molecule that is on top
proc plumedVis::findMolecule {args} {
   # Get the new molid
   foreach m [molinfo list] {
      if { [molinfo $m get top ] } { set topmol $m ; molinfo $m set drawn 1 }
      # There should only ever be one molecule drawn at a time (otherwise we can't make sense of colvars)
      if { [molinfo $m get drawn] && ![molinfo $m get top] } { molinfo $m set drawn 0 } 
   }
   return $topmol
}

# These routines read various types of file
proc plumedVis::browse {args} {
   variable filename

   set title [getargs $args "-title" {}]
   set todo [getargs $args "-todo" {}]

   set tmp [ tk_getOpenFile -initialdir [pwd] -title $title ]

   if { $todo=="readcolvars" && $tmp != "" } {
      set filename $tmp 
      readColvar 
   } elseif { $todo=="calccolvars" && $tmp != "" } {
      set filename $tmp
      openLenUnitDialogue
   } elseif { $todo=="readfes" && $tmp != "" } {
      set filename $tmp
      readFes -kT 0
   } elseif { $todo=="findbasins" && $tmp != "" } {
      set filename $tmp
      openBasinsDialogue
   } elseif { $todo=="sumhills" && $tmp != "" } {
      set filename $tmp
      openHillsDialogue
   } elseif { $todo=="readonions" && $tmp != "" } {
      set filename $tmp
       readOnions
   } else {
      puts "No mode specified for after in browse"
      return
   }
}

proc plumedVis::readOnions {args} {
   variable filename 
   variable cvdata
   variable recondata
   variable molid

   # Check for reconnaissance data
   set flag 0
   foreach cvname $cvdata([list $molid cv_names]) {
      set flag [expr $flag + $recondata([ list $molid $cvname isdata]) ] 
   }

   # Onions are only read when there is basin information already
   if { $flag==0 } {
      tk_messageBox -icon error -type ok -title Message -message "No basin information available for this molecule"
      return
   }

   set basmol $recondata([list $molid basin_molid])

   set fd [open $filename r]     ; # Open the ONIONS file 
   while { [gets $fd line] != -1 } {
      set sline [regexp -inline -all -- {\S+} $line]
      set bas [lindex $sline 1]    ; # the basin number is the second element on the line
      set cursize [lindex $cvdata([list $basmol $bas sizes]) end]
      # Read the size from the input file
      set sz [lindex $sline 5]     ; # the size is the sixth elelment on the line
      # Here the assumption is that the point will most likely be equidistant (in the metric) in all directions from the center
      set size [expr sqrt( $sz*$sz / $recondata([list $molid ncv]) )]   
      if { $size>$cursize } { 
         lappend cvdata([list $basmol $bas sizes]) $size
      } elseif { $size<[expr $cursize - 0.1] } {
         puts "ERROR - something has gone badly wrong $size $cursize"
      }
   }
   close $fd  
}

proc plumedVis::openLenUnitDialogue {args} {
   variable lenunit
   variable tfilename
    
   set tfilename "none"

   # Open len dialogue
   set len [toplevel ".len"]
   wm title $len "length units"
   wm resizable $len 0 0

   frame $len.top -padx 1m -pady 1m
   pack [label $len.lab -text "Please specify the length units"] -in $len.top -side top -fill both 
   frame $len.up -padx 1m -pady 1m
   label $len.up.lab -text "units" 
   menubutton $len.up.m -relief raised -bd 2 -direction flush -width 5 \
              -textvariable plumedVis::lenunit -menu $len.up.m.menu
   menu $len.up.m.menu
   grid $len.up.lab -row 1 -column 1 -sticky e
   grid $len.up.m   -row 1 -column 2 -sticky e
   pack $len.up -in $len.top -side top -fill both

   pack [label $len.flab -text "topology file containing masses and charges"] -in $len.top -side top -fill both
   frame $len.file -padx 1m -pady 1m
   entry $len.file.fent -width 5 -textvariable plumedVis::tfilename
   button $len.file.browse -text "Browse" -relief raised -command {set plumedVis::tfilename [tk_getOpenFile -initialdir [pwd] -title "Select the topology file"] } 
   grid $len.file.fent -row 1 -column 1 -sticky e
   grid $len.file.browse -row 1 -column 2 -sticky e
   pack $len.file -in $len.top -side top -fill both 

   # Setup the menu
   $len.up.m.menu delete 0 end
   $len.up.m configure -state disabled
   $len.up.m configure -state normal
   $len.up.m.menu add radiobutton -label "Angstroms" -value 1.0 -variable plumedVis::lenunit
   $len.up.m configure -state normal
   $len.up.m.menu add radiobutton -label "nm" -value 10 -variable plumedVis::lenunit
   $len.up.m configure -state normal
   $len.up.m.menu add radiobutton -label "au" -value 0.5292 -variable plumedVis::lenunit

   # Setup the buttons
   frame $len.but -padx 3m -pady 3m
   pack [button $len.but.cancel -text "Cancel" -relief raised -command {destroy .len ; return} ] -in $len.but -side right
   pack [button $len.but.ok -text "OK" -relief raised -command {[namespace code plumedVis::readPlumed] ; destroy .len } ] -in $len.but -side right
   pack $len.but -in $len.top -side top -fill both
   pack $len.top -fill both

   # Assume that lenghts are in Angstroms
   set lenunit 1.0
}

proc plumedVis::openHillsDialogue {args} {
   variable sumhills
   variable filename
   global env

   # Check sum hills is present  
   if { ![file exists "$env(plumedir)/utilities/sum_hills/sum_hills"] } {
      tk_messageBox -icon error -type ok -title Message -message "Could not find sum_hills should be in $env(plumedir)/utilities/sum_hills/sum_hills"   
      return
   }

   set fd [open $filename r]                        ; # Open the HILLS file
   gets $fd line                                    ; # Read the first line
   close $fd                                        ; # Close the HILLS file
   set sline [regexp -inline -all -- {\S+} $line]   ; # Convert the line to a list
   set nline [llength $sline]                       ; # Count the number of elements
   set ncv [expr ( $nline - 3 ) / 2 ]               ; # Calculate the number of cvs  

   # Check for integral ncv and that there enough cvs in the input COLVAR
   if { [expr $ncv - int($ncv)] > 0 } {
      tk_messageBox -icon error -type ok -title Message -message "Found a non-integral number of cvs in HILLS file abandoning"  
      return
   }
   set sumhills(ncv) $ncv        ; # Store number of cvs
   set sumhills(stride) 1000     ; # A tentative stride
   set sumhills(T) 300           ; # A tentative temperature 
   set sumhills(units) 0         ; # Stupid value for units

   # Open hills dialogue
   set hills [toplevel ".hills"]
   wm title $hills "Sum Hills"
   wm resizable $hills 0 0

   frame $hills.top -padx 1m -pady 1m

   frame $hills.fmess -padx 5m -pady 5m
   pack [label $hills.mess -text "Please specify the all information required to sum hills" -anchor s] -in $hills.fmess -side top
   pack $hills.fmess -in $hills.top -side top -fill both
 
   # Stride and kT
   frame $hills.mid -padx 1m -pady 1m
   labelframe $hills.dat -relief ridge -bd 2 -text "General" -padx 2m -pady 2m
   label $hills.dat.lstride -text "stride" -anchor e
   entry $hills.dat.stride -width 5 -textvariable plumedVis::sumhills(stride)
   label $hills.dat.ltemp -text "temperature" -anchor e
   entry $hills.dat.temp -width 5 -textvariable plumedVis::sumhills(T)
   label $hills.dat.lunits -text "units" -anchor e
   menubutton $hills.dat.units -relief raised -bd 2 -direction flush -width 5 \
             -textvariable plumedVis::sumhils(units) -menu $hills.dat.units.menu
   menu $hills.dat.units.menu
   grid $hills.dat.lstride -row 1 -column 1 -sticky e
   grid $hills.dat.stride  -row 1 -column 2 -sticky e
   grid $hills.dat.ltemp   -row 2 -column 1 -sticky e
   grid $hills.dat.temp    -row 2 -column 2 -sticky e
   grid $hills.dat.lunits  -row 3 -column 1 -sticky e
   grid $hills.dat.units   -row 3 -column 2 -sticky e
   pack $hills.dat -in $hills.mid -side left -fill both

   # CV Stuff
   labelframe $hills.cvs -relief ridge -bd 2 -text "Select CVs" -padx 2m -pady 2m

   label $hills.cvs.title -text "cv number" -anchor e
   label $hills.cvs.activ -text "active" -anchor e 
   label $hills.cvs.period -text "period" -anchor e
   grid $hills.cvs.title -row 1 -column 1 -sticky e
   grid $hills.cvs.activ -row 1 -column 2 -sticky e
   grid $hills.cvs.period -row 1 -column 3 -sticky e
   for { set i 1 } { $i<=$ncv } { incr i } { 
      label $hills.cvs.t$i -text "cv$i" -anchor e
      set plumedVis::sumhills([list active $i]) 0
      checkbutton $hills.cvs.a$i -variable plumedVis::sumhills([list active $i])  
      menubutton $hills.cvs.p$i -relief raised -bd 2 -direction flush -width 5 \
             -textvariable plumedVis::sumhils([list period $i]) -menu $hills.cvs.p$i.menu
      menu $hills.cvs.p$i.menu
      $hills.cvs.p$i.menu delete 0 end
      $hills.cvs.p$i configure -state disabled 
      $hills.cvs.p$i configure -state normal 
      $hills.cvs.p$i.menu add radiobutton -label "none" -value "none" -variable plumedVis::sumhills([list period $i])
      $hills.cvs.p$i configure -state normal
      $hills.cvs.p$i.menu add radiobutton -label "pi" -value pi -variable plumedVis::sumhills([list period $i]) 
      $hills.cvs.p$i configure -state normal
      $hills.cvs.p$i.menu add radiobutton -label "2pi" -value 2pi -variable plumedVis::sumhills([list period $i])
      grid $hills.cvs.t$i -row [expr $i+1] -column 1 -sticky e
      grid $hills.cvs.a$i -row [expr $i+1] -column 2 -sticky e
      grid $hills.cvs.p$i -row [expr $i+1] -column 3 -sticky e
   }
   pack $hills.cvs -in $hills.mid -side right -fill both
   pack $hills.mid -in $hills.top -side top -fill both

   # Setup the buttons
   frame $hills.but -padx 3m -pady 3m
   pack [button $hills.but.cancel -text "Cancel" -relief raised -command {destroy .hills ; return} ] -in $hills.but -side right
   pack [button $hills.but.ok -text "Sum Hills" -relief raised -command {[namespace code plumedVis::sumhills] ; destroy .hills } ] -in $hills.but -side right
   pack $hills.but -in $hills.top -side top -fill both
   pack $hills.top -fill both

   # Setup kT menu
   $hills.dat.units.menu delete 0 end
   $hills.dat.units configure -state disabled
   $hills.dat.units configure -state normal
   $hills.dat.units.menu add radiobutton -label "kJ/mol" -value 0.00831447 -variable plumedVis::sumhills(units)
   $hills.dat.units configure -state normal
   $hills.dat.units.menu add radiobutton -label "kcal/mol" -value 0.001987191 -variable plumedVis::sumhills(units)
   $hills.dat.units configure -state normal
   $hills.dat.units.menu add radiobutton -label "eV" -value 0.00008617343 -variable plumedVis::sumhills(units)
   $hills.dat.units configure -state normal
   $hills.dat.units.menu add radiobutton -label "dlpoly" -value 0.831451115 -variable plumedVis::sumhills(units)
   $hills.dat.units configure -state normal
   $hills.dat.units.menu add radiobutton -label "Rydbergs" -value 0.0000063363125 -variable plumedVis::sumhills(units)
   $hills.dat.units configure -state normal
   $hills.dat.units.menu add radiobutton -label "Natural" -value 1.0 -variable plumedVis::sumhills(units)

   # Set all periods to none by default
   for { set i 1 } { $i<=$ncv } { incr i } { set sumhills([list period $i]) "none" }
}

proc plumedVis::sumhills {args} {
   variable sumhills
   variable filename
   global env

   # Setup some stuff for calculations
   set kT [expr $sumhills(T)*$sumhills(units)]
   set period {}   ; # Set period to an empty list
   for { set i 1 } { $i<=$sumhills(ncv) } { incr i } {
      if { $sumhills([list active $i])==1 } { lappend ndw $i }
      if { $sumhills([list period $i])=="2pi" } {
         lappend period "-2pi $i"
      } elseif {$sumhills([list period $i])=="pi" } {
         lappend period "-pi $i"
      }
   }

   # These are all our sanity checks
   if { $kT==0 } {
      tk_messageBox -icon error -type ok -title Message -message "kT was not set correctly - check units" 
      return
   } elseif { ![string is integer -strict $sumhills(stride)] } {
      tk_messageBox -icon error -type ok -title Message -message "stride is not an integer - please remedy" 
      return
   } elseif { [llength $ndw] != 2 && [llength $ndw] != 1 } {
      tk_messageBox -icon error -type ok -title Message -message "Please select one or two cvs only" 
      return
   }

   # Setup a tempory directory to store fes.dat files
   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]" 

   # Copy HILLS file
   file copy -force $filename $tmpd
   set owd [pwd]
   cd $tmpd

   # Get the name of the hills file (i.e. without the path)
   set namelist [ lassign [split $filename /] ]
   set pos [ expr [llength $namelist] - 1 ] 
   set fname [lindex $namelist $pos]

   puts "Executing: $env(plumedir)/utilities/sum_hills/sum_hills -ndim $sumhills(ncv) -ndw [join $ndw] [join $period] -stride $sumhills(stride) -kt $kT -file $fname -out fes.dat" 
   catch {eval exec $env(plumedir)/utilities/sum_hills/sum_hills -ndim $sumhills(ncv) -ndw [join $ndw] [join $period] -stride $sumhills(stride) -kt $kT -file $fname -out fes.dat } sumhills_stdout

   # Check it worked
   if { ![file exists $tmpd/fes.dat] } {
      puts $sumhills_stdout
      puts "-------------"
      tk_messageBox -icon error -type ok -title Message -message "Something went wrong when running sum_hills.  Please see messages in log"
      return
   }

   # Read the fes
   set filename "$tmpd/fes.dat"   ;   readFes -kT $kT
   cd $owd   ; # Change back to current directory
}

proc plumedVis::openBasinsDialogue {args} {
   variable readcolvars
   variable molid
   variable cvdata
   variable recondata
   variable filename

   # First check that there is some colvar information available
   if { $readcolvars==0 } {
      tk_messageBox -icon error -type ok -title Message -message "No colvar information available - cannot find basins"    
      return
   }

   set fd [open $filename r]                        ; # Open the BASINS file 
   gets $fd line                                    ; # Read the first line
   close $fd                                        ; # Close the BASINS file
   set sline [regexp -inline -all -- {\S+} $line]   ; # Convert the line to a list
   set nline [llength $sline]                       ; # Count the number of elements
   set ncv [expr ( sqrt(4*$nline - 11) - 1 ) / 2 ]   ; # Calculate the number of cvs  

   # Check for integral ncv and that there enough cvs in the input COLVAR
   if { [expr $ncv - int($ncv)] > 0 } {
      tk_messageBox -icon error -type ok -title Message -message "Found a non-integral number of cvs in basins file abandoning"  
      return
   } elseif { $ncv > [llength [lsearch -all $cvdata([ list $molid cv_names ]) cv* ] ] } {
      tk_messageBox -icon error -type ok -title Message -message "Not enough calculated cvs to match  number of cvs in BASINS file"  
      return
   } 

   # Save the number of reconnaissance collective coordinates
   set recondata([list $molid ncv]) $ncv
   # Get the number of colvars we have data on
   set ncolvar [llength [lsearch -all $cvdata([ list $molid cv_names ]) cv* ] ]

   # Set the number of columns in cv matrix
   set ncols 3
   # Calculate the number of rows
   set nrows [expr int( $ncolvar / $ncols ) + 1]
   
   # Store the names of all the cvs in an array
   set n 0
   foreach cvname [lsearch -all -inline $cvdata([ list $molid cv_names ]) cv* ] {
      set cv_names($n) $cvname
      incr n
   }
   set nfcv $n   ; # Keep track of the total number of cvs

   # Open basins dialogue
   set bas [toplevel ".bas"]
   wm title $bas "Find Basins"
   wm resizable $bas 0 0

   frame $bas.top -padx 1m -pady 1m

   frame $bas.fmess -padx 5m -pady 5m
   pack [label $bas.mess -text "Please specify the $ncv collective coordinates used in reconnaissance" -anchor s] -in $bas.fmess -side top
   pack $bas.fmess -in $bas.top -side top -fill both

   # This creates all our little cv checkboxes - moment of gloating here this is so cool
   set n 0 
   frame $bas.fcv -padx 3m -pady 3m  
   for { set i 0 } { $i<$nrows } { incr i } {
      for { set j 0 } { $j<$ncols } { incr j } {
          set n [expr $i*$ncols + $j]
          if { $n>=$ncolvar } { break }
          label $bas.fcv.l$cv_names($n) -text "$cv_names($n)"
          # If the number of reconnaissance cvs is equal to the number of cvs there is data pre-select all cvs.
          if { $nfcv== $ncv } { set recondata([list $molid $cv_names($n) isdata]) 1 }
          checkbutton $bas.fcv.b$cv_names($n) -variable plumedVis::recondata([list $molid $cv_names($n) isdata])
          label $bas.fcv.sl$cv_names($n) -text "scale"
          set recondata([list $molid $cv_names($n) scalep]) 1.0
          # Scalep here allows us to scale variable where units in BASINS and vmd are different
          entry $bas.fcv.s$cv_names($n) -width 5 -textvariable plumedVis::recondata([list $molid $cv_names($n) scalep])
          grid $bas.fcv.l$cv_names($n) -row [expr $i+1] -column [expr 4*$j+1]
          grid $bas.fcv.b$cv_names($n) -row [expr $i+1] -column [expr 4*$j+2]
          grid $bas.fcv.sl$cv_names($n) -row [expr $i+1] -column [expr 4*$j+3]
          grid $bas.fcv.s$cv_names($n) -row [expr $i+1] -column [expr 4*$j+4]

      }
      if { $n>=$ncolvar } { break }
   }         

   # And the buttons that make it go
   frame $bas.but -padx 3m -pady 3m
   pack [button $bas.but.cancel -text "Cancel" -relief raised -command {destroy .bas ; return} ] -in $bas.but -side right
   pack [button $bas.but.ok -text "Find basins" -relief raised -command {destroy .bas ; [namespace code plumedVis::findBasins]} ] -in $bas.but -side right

   pack $bas.fcv -in $bas.top -side top -fill both
   pack $bas.but -in $bas.top -side top -fill both
   pack $bas.top -fill both
}

proc plumedVis::findBasins {args} {
   variable molid
   variable cvdata
   variable recondata
   variable recon_init_size
   variable filename

   global vmd_frame

   set n 0
   foreach cvname [lsearch -all -inline $cvdata([ list $molid cv_names ]) cv* ] {
      set n [expr $n + $recondata([list $molid $cvname isdata])]
   }

   if { $n != $recondata([list $molid ncv]) } {
      tk_messageBox -icon error -type ok -title Message -message "Incorrect number of collective coordinates specified"   
      openBasinsDialogue    ; # GAT test this doesn't do crazy stuff
      return
   }

   # We should now read in the basins
   set fd [open $filename r]     ; # Open the BASINS file
   set nbas 0                    ; # Intialize counter over number of basins

   # Now read in each line of the BASINS file
   while { [gets $fd line] != -1 } {
     set sline [regexp -inline -all -- {\S+} $line] 
     for {set j 0 } { $j<$recondata([list $molid ncv]) } { incr j } {
        lappend basins($nbas) [ lindex $sline [expr $j+2] ]
     }
     # Store the covariance data
     for {set j 0 } { $j<$recondata([list $molid ncv]) } { incr j } {
         for {set k 0 } { $k<$recondata([list $molid ncv]) } { incr k } {
             set ind [expr int( 2 + ( $j + 1 ) * $recondata([list $molid ncv]) + $k ) ]
             set basins([list $nbas $j $k]) [ lindex $sline $ind ]
         }
     }
     incr nbas
   }
   close $fd                                    ; # Close the BASINS file

   # This is the loop to find the nearest frame
   set nfr [molinfo $molid get numframes]
   for {set i 0} { $i<$nfr } { incr i } {

       # Calculate the distance from every basin
       for {set j 0} { $j<$nbas } { incr j } {
          if { $i==0 } { set mindist($j) 0 }
          set dist 0 ; set k 0
          foreach cvname [lsearch -all -inline $cvdata([ list $molid cv_names ]) cv* ] {
             if { $recondata([list $molid $cvname isdata])==1 } {
                 set tmp  [ expr $recondata([list $molid $cvname scalep])*[lindex $basins($j) $k] - [lindex $cvdata([list $molid $cvname]) $i] ]  
                 set dist [expr $dist+$tmp*$tmp] 
                 incr k
             }
          }
          # Find the closest frame
          if { $i==0 || $dist<$mindist($j) } { set basfr($j) $i  ; set mindist($j) $dist }
       }
   } 

   # Change to the tempory directory
   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"
   set owd [pwd]
   cd $tmpd

   # Copy the frames
   for {set j 0} { $j<$nbas } { incr j } {
      #puts "Nearest frame to basin $j is frame $basfr($j) distance $mindist($j)"
      # Output coordinate information on this basin to a pdb file
      [ atomselect $molid all frame $basfr($j) ] writepdb basin_$j.pdb
   }

   # Read in basin information
   if { $nbas > 0 } {
      trace remove variable vmd_frame write [namespace code updateHighlight]

      # tell vmd to read in the data and load it as a new molecule (N.B. the selection will be the new top molecule)
      mol new basin_0.pdb type pdb
      for { set i 1 } { $i<$nbas } { incr i } { mol addfile basin_$i.pdb type pdb }

      set newmol [findMolecule]
      mol delrep 0 $newmol        ; # Delete the default representation of the new molecule (we copy those of the old)

      # Get the list of representations for the old molecule 
      set nreps [molinfo $molid get numreps]
      for { set i 0 } { $i<$nreps } { incr i } {
          set toget "{ rep $i} {selection $i} {color $i}"
          set rep [ molinfo $molid get $toget ]
          mol addrep $newmol 
          mol modstyle $i $newmol [lindex $rep 0]
          mol modselect $i $newmol [lindex $rep 1]
          mol modcolor $i $newmol [lindex $rep 2]
      }

      # Copy recon data to new molecule
      set recondata([list $newmol ncv]) $recondata([list $molid ncv])

      # Create cvdata for basins
      set k 0
      foreach cvname [lsearch -all -inline $cvdata([ list $molid cv_names ]) cv* ] {
         if { $recondata([list $molid $cvname isdata])==1 } {
            # Get names of cvs for basins
            lappend cvdata([ list $newmol cv_names ]) $cvname
            # There is basin information available for the new molecule
            set recondata([list $newmol $cvname isdata]) 1
            # store cv values in recondata 
            for {set j 0} { $j<$nbas } { incr j } {
               lappend cvdata([ list $newmol $cvname ]) [lindex $basins($j) $k]
            }
            # Store shape parameters for basins molecule
            set n 0
            foreach cvname2 [lsearch -all -inline $cvdata([ list $molid cv_names ]) cv* ] {
               if { $recondata([list $molid $cvname2 isdata])==1 } {
                  for {set j 0} { $j<$nbas } { incr j } {
                      lappend cvdata([ list shape $newmol $cvname $cvname2]) $basins([list $j $k $n])
                  }
                  incr n
               }
            }
            incr k
         }
      }
      # Store initial sizes for all basins
      set init_size_tmp [expr sqrt( $recondata([list $molid ncv]) - 1 ) + $recon_init_size]
      set init_size [expr sqrt( $init_size_tmp*$init_size_tmp / $recondata([list $molid ncv]) ) ]
      for {set j 0} { $j<$nbas } { incr j } {
          lappend cvdata([ list $newmol $j sizes ]) $init_size
      }

      set recondata([list $molid basin_molid]) $newmol
      set recondata([list $newmol basin_molid]) $newmol
      set molid $newmol    ;# The top molecule has changed on read in
      updateMenus          ;# Change the menus to reflect change in molid

      # Turn the trace back on
      trace add variable vmd_frame write [namespace code updateHighlight]

      # Clean up all the pdb files we have created  
      eval file delete [glob -dir [pwd] basin_*]
   }
   cd $owd
}

proc plumedVis::readPlumed {args} {
   variable w
   variable filename
   variable tfilename
   variable lenunit

   global env

   # Check that driver is where it is supposed to be 
   if { ![file exists "$env(plumedir)/utilities/driver/driver"] } {
      tk_messageBox -icon error -type ok -title Message -message "Could not find driver should be in $env(plumedir)/utilities/driver/driver"
      return
   } 

   # open a file and setup the tmpdir to which we are copying the file
   set fd [open $filename r]
   set tmpd "[ plumedVis::tmpdir ]/plumedvis.[pid]"
   set od [open "$tmpd/metad.dat" w]

   # We must print colvars at every trajectory step
   puts $od "PRINT W_STRIDE 1"

   while { [gets $fd line] != -1 } {
      set sline [regexp -inline -all -- {\S+} $line]

      # Deal with CV units and extra files
      if { !([lsearch $sline "COORD"] < 0) || !([lsearch $sline "HBONDS"] < 0) || \
           !([lsearch $sline "WATERBRIDGE"] < 0) || !([lsearch $sline "ELSTPOT"] < 0) } {
         set flag 0
         set newcomm {}
         foreach elem $sline {
            if { $flag==1 } {
               lappend newcomm [expr $elem * $lenunit]
               set flag 0
            } else {
               lappend newcomm $elem
               if { $elem=="R_0" || $elem=="D_0" } { set flag 1 }
            }
         }
         puts $od [join $newcomm]
      } elseif { ![lsearch $sline "ENERGY"] < 0 } {
         tk_messageBox -icon error -type ok -title Message -message "It is not possible to calculate the energy using driver"    
         return
      } elseif { !([lsearch $sline "ALPHARMSD"] < 0) || !([lsearch $sline "ANTIBETARMSD"] < 0) || !([lsearch $sline "PARABETARMSD"] < 0) } {
          set flag 0
          set newcomm {}
          foreach elem $sline {
            if { $flag==2 } {
               lappend newcomm [expr $elem * $lenunit]
               set flag 0
            } elseif { $flag==1 } {
               set flag 0
            } else {
               if { $elem=="ANGSTROM SCALE" } {
                  set flag 1
               } elseif { $elem=="R_0" } {
                  set flag 2
                  lappend newcomm $ele
               } else {
                  lappend newcomm $ele
               }
            }
          }
          puts $od [join $newcomm]
      } elseif { !([lsearch $sline "S_PATH"] < 0) || !([lsearch $sline "Z_PATH"] < 0) || !([lsearch $sline "TARGETED"] < 0) } {
          if { !([lsearch -all $sline "CMAP"] < 0 } {
             if { $lenunit==1.0 } {
                set pos [lsearch $sline "INDEX"]
                set fname [join [ lindex $sline [expr $pos+1 ] ] ]
                file copy -force $fname* $tmpd
                set pos [lsearch $sline "MAP"]
                set fname [join [ lindex $sline [expr $pos+1 ] ] ]
                file copy -force $fname* $tmpd
                set pos [lsearch $sline "GROUP"]
                if { !($pos < 0) } {
                   set fname [join [ lindex $sline [expr $pos+1 ] ] ]
                   file copy -force $fname* $tmpd
                }
             } else {
                tk_messageBox -icon error -type ok -title Message -message "Cannot interpret PATH CV with CMAP if plumed.dat file is not in Angstroms"    
                return
             }
          } else {
             set pos [lsearch $sline "FRAMESET"]
             set fname [join [ lindex $sline [expr $pos+1 ] ] ]
             file copy -force $fname* $tmpd
          }
          puts $od $line
      } elseif { !([lsearch $sline "BESPOKE"] < 0) } {
         # GAT note the problems that might occur with distance cvs and this if
         #     we are not careful and ensure to take this into account when we
         #     develop this further
         set pos [ lsearch $sline "FILENAME" ]
         if { $pos>0 } {
            set fname [join [ lindex $sline [expr $pos+1 ] ] ]
            file copy -force $fname $tmpd
         } else {
            tk_messageBox -icon error -type ok -title Message -message "Syntax for bespoke cv is wrong - no filename"    
            return
         }
         puts $od $line
      } elseif { !([lsearch $sline "PCA"] < 0) } {
         set pos [lsearch  $sline "FRAME" ]
         if { $pos>0 } {
            set fname [join [ lindex $sline [expr $pos+1 ] ] ]
            file copy -force $fname $tmpd
         } else {
            tk_messageBox -icon error -type ok -title Message -message "Syntax for pca cv is wrong - no frame"
            return  
         }
         set pos [lsearch  $sline "EIGENVEC" ]
         if { $pos>0 } {
            set fname [join [ lindex $sline [expr $pos+1 ] ] ]
            file copy -force $fname $tmpd
         } else {
            tk_messageBox -icon error -type ok -title Message -message "Syntax for pca cv is wrong - no eigenvec"
            return
         }
      } elseif { !([lsearch $sline "RDF"] < 0) } {
         set flag 0
         set newcomm {}
         set bounds {}
         foreach elem $sline {
            if { $flag!=0 } {
               lappend newcomm [expr $elem * $lenunit]
               if { $flag==2 } {
                  lappend bounds [expr $elem * $lenunit]
               } elseif { $flag==1 && [ llength $bounds ]==1 } {
                  lappend bounds [expr $elem * $lenunit]
               }
               set flag [expr $flag-1]
            } else {
               lappend newcomm $elem
               if { $elem=="RANGE" } {
                  set flag 2
               } elseif { $elem=="WIDTH" } {
                  set flag 1
               }
            }
         }
         puts $od [join $newcomm] 
      } elseif { [lsearch $sline "PRINT"]<0 && [lsearch $sline "RECONNAISSANCE"]<0 && \
           [lsearch $sline "ONIONS"]<0 && [lsearch $sline "BASINS"]<0 && \
           [lsearch $sline "CLUSTER"]<0 && [lsearch $sline "HILLS"]<0 && \
           [lsearch $sline "WELLTEMPERED"]<0 } {
         # CVS dealt with here  : RGYR    DISTANCE  TORSION    POSITION    MINDIST    ANGLE
         #                        DIPOLE  DIHCOR    ALPHABETA  RMSDTOR     PUCKERING  HELIX
         puts $od $line
         # GAT note that syntax for walls etc is not OK  -- this should be dealt with perhaps?
      }
      # GAT missing CVs   ::  CMAP
   }
   close $fd   ; close $od    ; # Close the files 

   # Copy the topology file if it is setup
   if { $tfilename!="none" && $tfilename!="" } { file copy -force $tfilename $tmpd/topol.pdb }

   # Change to the tmp directory and output the trajectory in dcd format for driver
   set owd [pwd]   ; cd $tmpd  ; animate write dcd temp.dcd waitfor all top

   # Write out a topology file if we haven't read one in   
   if { ![file exists topol.pdb] } { 
      tk_messageBox -icon error -type ok -title Message -message "Warning - you did not specify a pdb file containing masses and charges \n please check masses and charges in $tmpd/topol.pdb"
      [atomselect top all] set occupancy [ [atomselect top all] get mass ]        
      [atomselect top all] set beta [ [atomselect top all] get charge ]
      [atomselect top all] writepdb topol.pdb 
   }

   # Get the periodic boundary conditions from vmd or by asking user
   set pbc [PeriodicBoundaryCondition] 
   if { [lsearch $pbc -cell]<0 && [lsearch $pbc -nopbc]<0 } {            
      tk_messageBox -icon warning -type ok -title Message -message "Sorry I couldn't find pbcs inside your input trajectory \n I have output all the file to run driver in $tmpd \n Please run driver manually and load resulting COLVAR file"
      cd $owd
      return 
   }

   # Now run driver  
   puts "Executing: $env(plumedir)/utilities/driver/driver -dcd temp.dcd -pdb topol.pdb -plumed metad.dat $pbc"   
   catch {eval exec $env(plumedir)/utilities/driver/driver -dcd temp.dcd -pdb topol.pdb -plumed metad.dat $pbc} driver_stdout
   
   # Check it worked
   if { ![file exists COLVAR] } {
      puts $driver_stdout
      puts "-------------"
      tk_messageBox -icon error -type ok -title Message -message "Something went wrong when running driver.  Please see messages in log"       
      cd $owd
      return
   } 

   # Read colvars
   set filename "$tmpd/COLVAR" ; readColvar 
   cd $owd    ;  # Put everything back to normal
}

proc plumedVis::readFes {args} {
   variable w
   variable appearance
   variable fesdata
   variable filename
   variable allcolvars

   # Remove any trace to the fesframe (just in case)
   trace remove variable plumedVis::fesdata(fesframe) write [namespace code updateFes]
   # Remove all colvar plots as we can't plot fes's with those
   ::gtPlot::remove -datatag allcolvars ;  set allcolvars 0 

   # Get kT from args
   set kT [getargs $args "-kT" 0]

   # Store all the fes data that is available
   findExtraFess

   set fd [open $filename r]    ; # Open fes file
   set np 0  ; # Counter over number of poitns in file
   # Get first line of file (is this a one or 2d fes)
   gets $fd line ; set sline [regexp -inline -all -- {\S+} $line]
  
   # Establish whether this is a one or 2d fes
   if { [llength $sline]==3 } {
      set fesdim 2 ; incr np
      lappend x [lindex $sline 0] ; lappend y [lindex $sline 1] ; lappend z [lindex $sline 2]
   } elseif { [llength $sline]==2 } {
      set fesdim 1 ; incr np
      # Note we set y here to an empty list just so that later thing will work
      lappend x [lindex $sline 0] ; lappend z [lindex $sline 1] ; set y {}
   } elseif { [llength $sline]!=0 } {
       tk_messageBox -icon error -type ok -title Message -message "Found a weird line in fes"
       return
   }

   while {[gets $fd line] != -1} {
      set sline [regexp -inline -all -- {\S+} $line]
      if { [llength $sline]==3 && $fesdim==2 } {
         lappend x [lindex $sline 0]
         lappend y [lindex $sline 1]
         lappend z [lindex $sline 2]
         incr np
      } elseif { [llength $sline]==2 && $fesdim==1 } {
         lappend x [lindex $sline 0] ; lappend z [lindex $sline 1]
         incr np
      } elseif { [llength $sline]!=0 } {
          tk_messageBox -icon error -type ok -title Message -message "Found a weird line in fes"    
          return
      }
   }
   close $fd

   # Save everthing - don't read these in so that we can overwrite old data
   set fesdata([ list fes np ]) $np
   set fesdata([ list fes x ]) $x
   set fesdata([ list fes y ]) $y
   set fesdata([ list fes z ]) $z 
   set fesdata([ list fes dim]) $fesdim

   # This opens a pop up dialogue so the user can tell gtplot the value of kT
   ::gtPlot::setkT -kT $kT

   # Now plot the data
   if { $fesdata([list fes dim])==1 } {
      # Pass the data to gtplot
      ::gtPlot::readData -datatag free_energy -npoints $np -xdata $fesdata([ list fes x ]) -ydata $fesdata([ list fes z ]) -isfes 1
      # Hopefully this draws a line through the 1d free energy surface
      ::gtPlot::plot2d -datatag free_energy -stride 1 -color black -drawline
   } else {
      # Pass the data to gtplot
      ::gtPlot::readData -datatag free_energy -npoints $np -xdata $fesdata([ list fes x ]) -ydata $fesdata([ list fes y ]) -zdata $fesdata([ list fes z]) -isfes 1
      # This plots the data
      ::gtPlot::plot2d -datatag free_energy -fescolors $appearance(fescol) -ncolors $appearance(nfescolors)
      # This sets up the zoom (only needed for 2d surfaces)
      setupZoom
   }

   # And add a trace to the fesframe
   trace add variable plumedVis::fesdata(fesframe) write [namespace code updateFes]
}

# This allows us to plot the free energy as a function of simulation time
proc plumedVis::updateFes {args} {
   variable w
   variable appearance
   variable fesdata
   variable filename

   # Plot new fes
   if { $fesdata(fesframe)<=$fesdata(nfes) } {
      set fd [open $fesdata(fileroot).$fesdata(fesframe) r]    ; # Open fes file

      set np 0
      while {[gets $fd line] != -1} {
         set sline [regexp -inline -all -- {\S+} $line]
         if { [llength $sline]==3 && $fesdata([ list fes dim ])==2 } {
            lappend x [lindex $sline 0]
            lappend y [lindex $sline 1]
            if { [lindex $x $np] != [lindex $fesdata([ list fes x ]) $np] || [lindex $y $np] != [lindex $fesdata([ list fes y ]) $np] } {
               tk_messageBox -icon error -type ok -title Message -message "Grids in fes files do not match so cannot compare" 
               return
            }
            lappend z [lindex $sline 2]
            incr np
         } elseif { [llength $sline]==2 && $fesdata([ list fes dim ])==1 } {
            lappend x [lindex $sline 0]
            if { [lindex $x $np] != [lindex $fesdata([ list fes x ]) $np] } {
              tk_messageBox -icon error -type ok -title Message -message "Grids in fes files do not match so cannot compare"
              return
            }
            lappend z [lindex $sline 1]
            incr np
         } elseif { [llength $sline]!=0 } {
             tk_messageBox -icon error -type ok -title Message -message "Found a weird line in fes"    
             return
         }
      }
      close $fd

      # Do we want to plot the difference from the final fes
      if { $appearance(plotdiff)==1 } {
         for {set i 0} {$i<$np} {incr i} {
            lappend pz [ expr [lindex $fesdata([ list fes z ]) $i] - [lindex $z $i] ]
         }
      } else {
         set pz $z
      }

      # Do we want to rescale the z axis
      if { $appearance(fixz)==1 } {
         set fixztag "-nonewbounds"
      } else {
         set fixztag {}
      }

      if { $fesdata([list fes dim])==1 } {
        # This passes 1D fes data to gtplot
        ::gtPlot::readData -datatag free_energy -npoints $np -xdata $x -ydata $pz -isfes 1
        # This plots the data
        ::gtPlot::plot2d -datatag free_energy -stride 1 -color black -drawline
      } else {
        # This reads the data from input
        ::gtPlot::readData -datatag free_energy -npoints $np -xdata $x -ydata $y -zdata $pz $fixztag -isfes 1
        # This plots the data
        ::gtPlot::plot2d -datatag free_energy -fescolors $appearance(fescol) -ncolors $appearance(nfescolors)
      }

   } elseif { $fesdata(fesframe)<1 } {
      tk_messageBox -icon error -type ok -title Message -message "There are no fes files with an index less than 1" 
      return 
   } else {
      set np $fesdata([list fes np])
      if { $fesdata([list fes dim])==1 } {
        # This passes 1D fes data to gtplot
        ::gtPlot::readData -datatag free_energy -npoints $np -xdata $fesdata([ list fes x ]) -ydata $fesdata([ list fes z ]) -isfes 1
        # This plots the data
        ::gtPlot::plot2d -datatag free_energy -stride 1 -color black -drawline
      } else {
         # This reads the data from input
         ::gtPlot::readData -datatag free_energy -npoints $np -xdata $fesdata([ list fes x ]) -ydata $fesdata([ list fes y ]) -zdata $fesdata([ list fes z]) -isfes 1
         # This plots the data
         ::gtPlot::plot2d -datatag free_energy -fescolors $appearance(fescol) -ncolors $appearance(nfescolors)
      }
   } 
}

proc plumedVis::findExtraFess {args} {
   variable w
   variable filename
   variable fesdata

   set fesdata(fileroot) $filename

   # This getst the filename only
   set namelist [ lassign [split $filename /] ]
   set pos [ expr [llength $namelist] - 1 ]
   set fname [lindex $namelist $pos]

   # This get the directory path
   for {set i 0} {$i<$pos} {incr i} {
      lappend tpath [lindex $namelist $i]
   } 
   set path [join $tpath /]

   set fesdata(nfes) [expr [llength [glob -dir $path $fname*] ] - 1] 

   # This sets up the spinbox menu so that we can scroll through frames
   set fesdata(frames) {}
   for {set i 1} {$i<=$fesdata(nfes)} {incr i} {
      lappend fesdata(frames) $i
   }
   lappend fesdata(frames) [expr $fesdata(nfes) + 1]
   $w.top.pmenus.spinb configure -values $plumedVis::fesdata(frames) 

   # We are plotting the final free energy
   set fesdata(fesframe) [expr $fesdata(nfes) + 1]
}

# This calculate the free energy in the highlighted region/s as a function
# of time
proc plumedVis::drawConvergencePlot {args} {
   variable fesdata

   for {set i 1} {$i<=$fesdata(nfes)} {incr i} {

       # Clear the lists that hold the data points
       set x {}   ; set y {}   ; set z {}

       # Open and read the fes.dat.$i file 
       set fd [open $fesdata(fileroot).$i r]    

       set np 0
       while {[gets $fd line] != -1} {
          set sline [regexp -inline -all -- {\S+} $line]
          if { [llength $sline]==3 && $fesdata([list fes dim])==2 } {
             lappend x [lindex $sline 0]
             lappend y [lindex $sline 1]
             if { [lindex $x $np] != [lindex $fesdata([ list fes x ]) $np] || [lindex $y $np] != [lindex $fesdata([ list fes y ]) $np] } {
               tk_messageBox -icon error -type ok -title Message -message "Grids in fes files do not match so cannot compare"
               return
             }
             lappend z [lindex $sline 2]
             incr np
          } elseif { [llength $sline]==2 && $fesdata([list fes dim])==1 } {
             lappend x [lindex $sline 0]
             if { [lindex $x $np] != [lindex $fesdata([ list fes x ]) $np] } {
               tk_messageBox -icon error -type ok -title Message -message "Grids in fes files do not match so cannot compare"
               return
             }
             lappend z [lindex $sline 1]
             incr np
          } elseif { [llength $sline]!=0 } {
              tk_messageBox -icon error -type ok -title Message -message "Found a weird line in fes"    
              return
          }
       }
       close $fd 

       # Pass the data to gtPlot so that it can work out free energy differences
       if { $fesdata([list fes dim])==1 } { 
          ::gtPlot::readData -datatag tfes -npoints $np -xdata $x -ydata $z
       } else {
          ::gtPlot::readData -datatag tfes -npoints $np -xdata $x -ydata $y -zdata $z
       }
       # Calculate the free energy defined by the highlight for this free energy surface
       lappend fesval [::gtPlot::calcFes -festag free_energy -datatag tfes]
       lappend t $i

   }
   # This is the difference for the final calculated fes
   lappend fesval [::gtPlot::calcFes -festag free_energy -datatag free_energy]
   lappend t [expr $fesdata(nfes)+1] 

   # And plot the data
   set convplot [multiplot -title "Free energy as a function of time" -xlabel "Frame" -nostats]

   $convplot add $t $fesval -legend "Free energy" -lines -marker circle \
        -radius 2 -fillcolor blue -color blue -nostats 
   $convplot replot
}

proc plumedVis::setupZoom {args} {
   variable appearance
   ::gtPlot::setproperty -datatag free_energy -fesZoomStyle $appearance(fixz)
}

proc plumedVis::replotFes {args} {
   variable appearance

   # This plots the data 
   ::gtPlot::plot2d -datatag free_energy -fescolors $appearance(fescol) -ncolors $appearance(nfescolors)
}

proc plumedVis::readColvar {args} {
   variable w
   variable appearance
   variable molid
   variable xcoord
   variable ycoord
   variable cvdata
   variable recondata
   variable filename
   variable readcolvars

   set fd [open $filename r]                 ; # Open COLVAR FILE
   set ncv [ expr 0 ]                        ; # Counter over Column number
   gets $fd header                           ; # Get the first line of the COLVAR file
   set sheader [regexp -inline -all -- {\S+} $header] 
   
   set ncv [llength $sheader]                ; # Get the number of fields in the COLVAR file

   if { [lindex $sheader 1] != "FIELDS" } { 
      tk_messageBox -icon error -type ok -title Message -message "Input file is not a COLVAR file - cannot read"   
      return
   }

   # Read the header and create the menus
   for { set id 2 } { $id < $ncv } { incr id } {
     set cv_name [lindex $sheader $id]
     # We must store the list of cv names so that we can copy colvar information 
     # when we select datapoints
     lappend cvdata([list $molid cv_names]) $cv_name
     # This tells us whether there is basin information for this colvar name
     # We set it up here but it only comes into its own once we load basins
     set recondata([ list $molid $cv_name isdata]) 0  
   }

   # Read the colvar file
   set ncv [expr $ncv - 2]    ; # Don't count the first two fields in the colvar header
   set ncolvar 0              ; # Counter of the number of colvars read in
   while {[gets $fd line] != -1} {
      set sline [regexp -inline -all -- {\S+} $line]
      for { set id 0 } { $id < $ncv } { incr id } {
          set cv_name [ lindex $cvdata([list $molid cv_names]) $id ]
          lappend cvdata([ list $molid $cv_name ]) [lindex $sline $id ]
      }
      incr ncolvar
   }
   close $fd

   if { $ncolvar != [molinfo $molid get numframes] } {
       tk_messageBox -icon error -type ok -title Message -message "Number of frames in colvar does not match number of frames in top trajectory - aborting"  
       return
   } else {
       updateMenus  ;           # Setup the cv menus
       set readcolvars 1 ;      # We have colvar information
  }
}

proc plumedVis::updateMenus {args} {
   variable molid
   variable cvdata
   variable xcoord
   variable ycoord
   variable w

   if { $molid != [findMolecule] } {
      puts "ERROR -- Something has gone wrong in updating a menu somewhere - fix me Gareth"
      return
   }

   # Clear old cv menus
   $w.top.pmenus.xaxis.m.menu delete 0 end
   $w.top.pmenus.yaxis.m.menu delete 0 end
   $w.top.pmenus.xaxis.m configure -state disabled
   $w.top.pmenus.yaxis.m configure -state disabled

   # Create none switches
   $w.top.pmenus.xaxis.m configure -state normal
   $w.top.pmenus.xaxis.m.menu add radiobutton -label "none" -value "none" -variable plumedVis::xcoord
   $w.top.pmenus.yaxis.m configure -state normal
   $w.top.pmenus.yaxis.m.menu add radiobutton -label "none" -value "none" -variable plumedVis::ycoord

   # Now create the menus
   foreach cvname $cvdata([list $molid cv_names]) {
      $w.top.pmenus.xaxis.m configure -state normal 
      $w.top.pmenus.xaxis.m.menu add radiobutton -value $cvname \
          -label "$cvname" -variable plumedVis::xcoord
      $w.top.pmenus.yaxis.m configure -state normal 
      $w.top.pmenus.yaxis.m.menu add radiobutton -value $cvname \
          -label "$cvname" -variable plumedVis::ycoord
   }

   # We currently have not selected xcoord or ycoord
   set xcoord "none"    ; set ycoord "none"
}

# These are the tools for drawing the parts of the window
proc plumedVis::createMenubar {w} {
  frame $w -relief raised -bd 2

  # File menu
  menubutton $w.file -text "File   " -underline 0 -menu $w.file.menu
  $w.file config -width 5
  pack $w.file -side left
  menu $w.file.menu -tearoff no
  $w.file.menu add command -label "Load colvars" -command [namespace code {plumedVis::browse -todo "readcolvars" -title "Please select the COLVAR file to read in"} ]
  $w.file.menu add command -label "Load fes" -command [namespace code {plumedVis::browse -todo "readfes" -title "Please select the 2d FES to read in"} ]
  $w.file.menu add command -label "Load onions" -command [namespace code {plumedVis::browse -todo "readonions" -title "Select the ONIONS file from a recon metad simuiation"} ]
  $w.file.menu add command -label "Print graph" -command [namespace code {::gtPlot::printCanvas $w} ]  ; # GAT this command needs testing   

  # Calculate menu
  menubutton $w.calc -text "Calculate   " -underline 0 -menu $w.calc.menu
  $w.calc config -width 10
  pack $w.calc -side left
  menu $w.calc.menu -tearoff no
  $w.calc.menu add command -label "Calculate colvars" -command [namespace code {plumedVis::browse -todo "calccolvars" -title "Please select a plumed input file"} ]
  $w.calc.menu add command -label "Save selection" -command [namespace code plumedVis::storeSelection]
  $w.calc.menu add command -label "Sum hills" -command [namespace code {plumedVis::browse -todo "sumhills" -title "Please select a HILLS file"} ]
  $w.calc.menu add command -label "Convergence of selected free energy" -command [namespace code plumedVis::drawConvergencePlot]
  $w.calc.menu add command -label "Find reconnaissance basins" -command [namespace code {plumedVis::browse -todo "findbasins" -title "Please select a BASINS file"} ]

  menubutton $w.look -text "Appearance   " -underline 0 -menu $w.look.menu
  $w.look config -width 15
  pack $w.look -side left
  menu $w.look.menu -tearoff no
  $w.look.menu add command -label "Colvar appearance" -command [namespace code plumedVis::openCVlook]
  $w.look.menu add command -label "Fes appearance" -command [namespace code plumedVis::openfeslook]

  # Help menu
  menubutton $w.help -text "Help   " -menu $w.help.menu
  $w.help config -width 5
  pack $w.help -side right
  menu $w.help.menu -tearoff no
  $w.help.menu add command -label "Plumed Help..." -command "vmd_open_url http://merlino.mi.infn.it/~plumed/PLUMED/Home.html"

  return $w
}

proc plumedVis::openfeslook {args} {
  variable fes
  variable appearance

  if [winfo exists .feslook] {
    wm deiconify $fes
    return
  }

  # Get information on what the fes scale currently looks like from gtplot
  set stuff [ ::gtPlot::getFesScaleAppearance ]
  set appearance(zticspacing) [lindex $stuff 0]
  set appearance(znmintics) [lindex $stuff 1]   
  set appearance(zmajticl) [lindex $stuff 2]
  set appearance(zminticl) [lindex $stuff 3]

  set fes [toplevel ".feslook"]
  wm title $fes "FES Appearance"
  wm resizable $fes 0 0
  bind $fes <Destroy> plumedVis::destroyFeslook

  frame $fes.top -padx 1m -pady 1m

  labelframe $fes.afram -relief ridge -bd 2 -text "FES appearance" -padx 2m -pady 2m
  label $fes.afram.nlab -text "n. colors" -anchor e
  entry $fes.afram.nent -width 5 -textvariable plumedVis::appearance(nfescolors)
  label $fes.afram.clab -text "colorscale" -anchor e
  menubutton $fes.afram.c -relief raised -bd 2 -direction flush -textvariable plumedVis::appearance(fescol) -menu $fes.afram.c.menu
  menu $fes.afram.c.menu
  $fes.afram.c configure -state disabled

  $fes.afram.c.menu delete 0 end
  $fes.afram.c configure -state normal
  $fes.afram.c.menu add radiobutton -label jet -value jet -variable plumedVis::appearance(fescol)
  $fes.afram.c configure -state normal
  $fes.afram.c.menu add radiobutton -label hsv -value hsv -variable plumedVis::appearance(fescol)
  $fes.afram.c configure -state normal
  $fes.afram.c.menu add radiobutton -label hot -value hot -variable plumedVis::appearance(fescol)
  $fes.afram.c configure -state normal
  $fes.afram.c.menu add radiobutton -label cool -value cool -variable plumedVis::appearance(fescol)
  $fes.afram.c configure -state normal
  $fes.afram.c.menu add radiobutton -label grey -value grey -variable plumedVis::appearance(fescol)

  grid $fes.afram.nlab -row 1 -column 1 -sticky e
  grid $fes.afram.nent -row 1 -column 2 -sticky e
  grid $fes.afram.clab -row 2 -column 1 -sticky e
  grid $fes.afram.c -row 2 -column 2 -sticky e
  pack $fes.afram -in $fes.top -side top -fill both

  labelframe $fes.sfram -relief ridge -bd 2 -text "Major scale tics" -padx 2m -pady 2m
  label $fes.sfram.spal -text "spacing" -anchor e
  entry $fes.sfram.spa -width 5 -textvariable plumedVis::appearance(zticspacing)
  label $fes.sfram.len -text "tic length" -anchor e
  eval spinbox $fes.sfram.en -textvariable plumedVis::appearance(zmajticl) -from 1.0 -to 8.0 -increment 0.25 -width 5
  grid $fes.sfram.spal -row 1 -column 1 -sticky e
  grid $fes.sfram.spa -row 1 -column 2 -sticky e
  grid $fes.sfram.len  -row 2 -column 1 -sticky e
  grid $fes.sfram.en   -row 2 -column 2 -sticky e
  pack $fes.sfram -in $fes.top -side top -fill both

  labelframe $fes.mfram -relief ridge -bd 1 -text "Minor scale Tics" -padx 1m -pady 1m
  label $fes.mfram.nl -text "n. tics" -anchor e
  eval spinbox $fes.mfram.n -textvariable plumedVis::appearance(znmintics) -from 0 -to 5 -increment 1 -width 5
  label $fes.mfram.lenl -text "tic length" -anchor e
  eval spinbox $fes.mfram.len -textvariable plumedVis::appearance(zminticl) -from 1.0 -to 4.0 -increment 0.125 -width 5
  grid $fes.mfram.nl   -row 1 -column 1 -sticky e  
  grid $fes.mfram.n    -row 1 -column 2 -sticky e
  grid $fes.mfram.lenl -row 2 -column 1 -sticky e 
  grid $fes.mfram.len  -row 2 -column 2 -sticky e
  pack $fes.mfram -in $fes.top -side top -fill both

  # Now create our buttons
  frame $fes.bfram -padx 3m -pady 3m
  pack [button $fes.bfram.ok -text "Replot" -relief raised -command {[namespace code plumedVis::replotFes]; [namespace code plumedVis::destroyFeslook]; destroy .feslook} ] -in $fes.bfram -side right
  pack [button $fes.bfram.cancel -text "Cancel" -relief raised -command {[namespace code plumedVis::destroyFeslook]; destroy .feslook} ] -in $fes.bfram -side right
  pack $fes.bfram -in $fes.top -side top -fill both
  pack $fes.top -fill both

  # We must set up some traces on the appearance of the fes scale
  trace add variable appearance(zticspacing) write [ namespace code plumedVis::redrawScale ]  
  trace add variable appearance(znmintics) write [ namespace code plumedVis::redrawScale ] 
  trace add variable appearance(zmajticl) write [ namespace code plumedVis::redrawScale ] 
  trace add variable appearance(zminticl) write [ namespace code plumedVis::redrawScale ] 
}

proc plumedVis::destroyFeslook {args} {
  variable appearance
  trace remove variable appearance(zticspacing) write [ namespace code plumedVis::redrawScale ]
  trace remove variable appearance(znmintics) write [ namespace code plumedVis::redrawScale ]  
  trace remove variable appearance(zmajticl) write [ namespace code plumedVis::redrawScale ]   
  trace remove variable appearance(zminticl) write [ namespace code plumedVis::redrawScale ]   
}

proc plumedVis::redrawScale {args} {
  variable appearance
  ::gtPlot::setFesScaleAppearance $appearance(zticspacing) $appearance(znmintics) $appearance(zmajticl) $appearance(zminticl)
}

proc plumedVis::openCVlook {args} {
  variable cv
  variable appearance
  variable highlight_nbefore
  variable highlight_nafter
  variable highlight_lines
  variable cvstride

  if [winfo exists .cvlook] {
    wm deiconify $cv
    return
  }

  set cv [toplevel ".cvlook"]
  wm title $cv "CV Plotter Appearance"
  wm resizable $cv 0 0
  bind $cv <Destroy> plumedVis::destroyCVlook

  frame $cv.top -padx 1m -pady 1m

  labelframe $cv.afram -relief ridge -bd 2 -text "plot appearance" -padx 2m -pady 2m
  label $cv.afram.lab -text "Plot stride" -anchor e 
  entry $cv.afram.ent -width 5 -textvariable plumedVis::cvstride 
  label $cv.afram.slab -text "Point size" -anchor e 
  eval spinbox $cv.afram.sspin -textvariable plumedVis::appearance(boxSize) -from 2 -to 30 -increment 1 -width 5 

  grid $cv.afram.lab -row 1 -column 1 -sticky e
  grid $cv.afram.ent -row 1 -column 2 -sticky e
  grid $cv.afram.slab -row 2 -column 1 -sticky e
  grid $cv.afram.sspin -row 2 -column 2 -sticky e
  pack $cv.afram -in $cv.top -side top -fill both

  labelframe $cv.hfram -relief ridge -bd 2 -text "highlight appearance" -padx 2m -pady 2m
  label $cv.hfram.lnb4 -text "Highlight frames before current"
  entry $cv.hfram.enb4 -width 5 -textvariable plumedVis::highlight_nbefore
  label $cv.hfram.lna4 -text "Highlight frames after current"
  entry $cv.hfram.ena4 -width 5 -textvariable plumedVis::highlight_nafter
  label $cv.hfram.clab -text "Colorscale for highlight"
  menubutton $cv.hfram.c -relief raised -bd 2 -direction flush -textvariable plumedVis::appearance(trajcol) -menu $cv.hfram.c.menu
  menu $cv.hfram.c.menu
  $cv.hfram.c configure -state disabled
  $cv.hfram.c.menu delete 0 end
  $cv.hfram.c configure -state normal
  $cv.hfram.c.menu add radiobutton -label jet -value jet -variable plumedVis::appearance(trajcol)
  $cv.hfram.c configure -state normal
  $cv.hfram.c.menu add radiobutton -label hsv -value hsv -variable plumedVis::appearance(trajcol)
  $cv.hfram.c configure -state normal
  $cv.hfram.c.menu add radiobutton -label hot -value hot -variable plumedVis::appearance(trajcol)
  $cv.hfram.c configure -state normal
  $cv.hfram.c.menu add radiobutton -label cool -value cool -variable plumedVis::appearance(trajcol)
  $cv.hfram.c configure -state normal
  $cv.hfram.c.menu add radiobutton -label grey -value grey -variable plumedVis::appearance(trajcol)

  label $cv.hfram.llla -text "Draw line through highlighted points"
  checkbutton $cv.hfram.lbox -variable plumedVis::highlight_lines

  grid $cv.hfram.lnb4 -row 1 -column 1 -sticky e
  grid $cv.hfram.enb4 -row 1 -column 2 -sticky w
  grid $cv.hfram.lna4 -row 2 -column 1 -sticky e
  grid $cv.hfram.ena4 -row 2 -column 2 -sticky w
  grid $cv.hfram.clab -row 3 -column 1 -sticky e
  grid $cv.hfram.c -row 3 -column 2 -sticky w
  grid $cv.hfram.llla -row 4 -column 1 -sticky e
  grid $cv.hfram.lbox -row 4 -column 2 -sticky w 
  pack $cv.hfram -in $cv.top -side top -fill both

  frame $cv.d -padx 3m -pady 3m
  pack [ button $cv.d.button  -text "Dismiss" -relief raised -command { [namespace code plumedVis::destroyCVlook]; destroy .cvlook} ] -in $cv.d -side right 
  pack $cv.d -in $cv.top -side top
  pack $cv.top -fill both 

  # These control the way the trajectory is displayed
  trace add variable plumedVis::appearance(boxSize) write [namespace code plotData]
  trace add variable plumedVis::cvstride write [namespace code plotData]
  trace add variable plumedVis::appearance(trajcol) write [namespace code updateHighlight]
  trace add variable plumedVis::highlight_nbefore write [namespace code updateHighlight]
  trace add variable plumedVis::highlight_nafter write [namespace code updateHighlight]
  trace add variable plumedVis::highlight_lines write [namespace code updateHighlight]

  return $cv
}

proc plumedVis::destroyCVlook {args} {
  variable appearance 
  variable cvstride
  variable highlight_nbefore 
  variable highlight_nafter 
  variable highlight_lines 

  trace remove variable plumedVis::appearance(boxSize) write [namespace code plotData]
  trace remove variable plumedVis::cvstride write [namespace code plotData]
  trace remove variable plumedVis::appearance(trajcol) write [namespace code updateHighlight]
  trace remove variable plumedVis::highlight_nbefore write [namespace code updateHighlight]
  trace remove variable plumedVis::highlight_nafter write [namespace code updateHighlight]
  trace remove variable plumedVis::highlight_lines write [namespace code updateHighlight]
}

proc plumedVis::labelspinbox {w args} {

   set text [getargs $args "-text" {}]
   # This removes the -text argument
   set pos [lsearch $args -text]
   if {$pos>=0} {
      set args [lreplace $args $pos [expr {$pos+1}]]
   }

   frame $w
   label $w.xlabel -text "$text "
   eval spinbox $w.dimx $args
   grid $w.xlabel -row 1 -column 1 -sticky we
   grid $w.dimx   -row 1 -column 2 -sticky we

   return $w
}

proc plumedVis::axisMenu {w args} {
   set text [getargs $args "-text" {}]
   set var  [getargs $args "-textvariable" {}]
   
   frame $w
   label $w.lab -text $text -anchor w
   menubutton $w.m -relief raised -bd 2 -direction flush -width 5 \
         -textvariable $var -menu $w.m.menu
   menu $w.m.menu  

   $w.m.menu delete 0 end
   $w.m configure -state disabled 

   # No cv to be loaded
   $w.m configure -state normal
   $w.m.menu add radiobutton -label "none" -value "none" -variable $var

   grid $w.lab -row 1 -column 1 -sticky we
   grid $w.m -row 1 -column 2 -sticky we

   return $w
}

# This gets arguments from a list
proc plumedVis::getargs {arglist tag defl {n 1}} {
   set pos [lsearch $arglist $tag]
   if {$pos<0}  { return $defl}
   return [join [ lrange $arglist [expr $pos+$n ] [expr $pos+$n ] ] ]
}

# This returns a tempory directory
proc plumedVis::tmpdir {args} {
   global tcl_platform
   switch $tcl_platform(platform) {
       unix {
           set tmpdir /tmp   ;# or even $::env(TMPDIR), at times.
       } macintosh {
           set tmpdir $::env(TRASH_FOLDER)  ;# a better place?
       } default {
           set tmpdir [pwd]
           catch {set tmpdir $::env(TMP)}
           catch {set tmpdir $::env(TEMP)}
       }
   }
}

# This returns pbc input for driver
proc plumedVis::PeriodicBoundaryCondition {args} {
  if { [molinfo top get alpha]!=90 || [molinfo top get beta]!=90 || [molinfo top get gamma]!=90 } {
     tk_messageBox -icon error -type ok -title Message -message "Can't run driver with non-orthorhombic unit cells"   
     return -1
  }
  if { [molinfo top get a]==0 } {
     set ans [tk_messageBox -icon question -type yesno -title Message -message "Are there really no pbcs?"]
     switch -- $ans {
         yes { return -nopbc }
         no  { return -1 }
     } 
  } else {
     return "-cell [molinfo top get a] [molinfo top get b] [molinfo top get c]"
  }
}
