###############################################################
#-- INITIALIZATION
#--
#-- Load Glyph package and prepare Pointwise. Also define
#-- the working directory, base input file name
#--
###############################################################

# Load Glyph
package require PWI_Glyph 2.18.0

# Setup Pointwise
pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}
pw::Application clearModified

# Working directory and base input CGNS file name
set cwd [file dirname [info script]]
set basename "hump2newtop_noplenumZ409x109"


###############################################################
#-- MAIN SCRIPT
#--
###############################################################

# -------------------------------------
# Read CGNS mesh
# -------------------------------------
set _TMP(mode_1) [pw::Application begin GridImport]
  $_TMP(mode_1) initialize -strict -type Automatic [file join $cwd $basename.cgns]
  $_TMP(mode_1) read
  $_TMP(mode_1) convert
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Import Grid}

# -------------------------------------
# Set the solver
# -------------------------------------
pw::Application setCAESolver {EXODUS II} 3
pw::Application markUndoLevel {Select Solver}

# -------------------------------------
# Set the boundary conditions
# -------------------------------------
set _DM(1) [pw::GridEntity getByName "dom-1"]
set _BL(1) [pw::GridEntity getByName "Zone   1"]
set _DM(2) [pw::GridEntity getByName "dom-2"]
set _DM(3) [pw::GridEntity getByName "dom-3"]
set _DM(4) [pw::GridEntity getByName "dom-4"]
set _DM(5) [pw::GridEntity getByName "dom-5"]
set _DM(6) [pw::GridEntity getByName "dom-6"]
set _TMP(PW_1) [pw::BoundaryCondition getByName "Unspecified"]
unset _TMP(PW_1)
set _TMP(PW_2) [pw::BoundaryCondition getByName "Unspecified"]
set _TMP(PW_3) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_4) [pw::BoundaryCondition getByName "bc-2"]
unset _TMP(PW_3)
set _TMP(PW_5) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_6) [pw::BoundaryCondition getByName "bc-3"]
unset _TMP(PW_5)
set _TMP(PW_7) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_8) [pw::BoundaryCondition getByName "bc-4"]
unset _TMP(PW_7)
set _TMP(PW_9) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_10) [pw::BoundaryCondition getByName "bc-5"]
unset _TMP(PW_9)
set _TMP(PW_11) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_12) [pw::BoundaryCondition getByName "bc-6"]
unset _TMP(PW_11)
$_TMP(PW_4) setName "bottomwall"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_6) setName "top"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_8) setName "inlet"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_10) setName "outlet"
pw::Application markUndoLevel {Name BC}

set _TMP(PW_13) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_14) [pw::BoundaryCondition getByName "bc-7"]
unset _TMP(PW_13)
$_TMP(PW_12) setName "front"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_14) setName "back"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_4) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_6) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_8) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_10) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_12) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_14) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_14) apply [list [list $_BL(1) $_DM(1)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_12) apply [list [list $_BL(1) $_DM(4)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_8) apply [list [list $_BL(1) $_DM(2)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_10) apply [list [list $_BL(1) $_DM(5)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_6) apply [list [list $_BL(1) $_DM(6)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_4) apply [list [list $_BL(1) $_DM(3)]]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_2)
unset _TMP(PW_4)
unset _TMP(PW_6)
unset _TMP(PW_8)
unset _TMP(PW_10)
unset _TMP(PW_12)
unset _TMP(PW_14)
set _TMP(PW_15) [pw::BoundaryCondition getByName "Unspecified"]
unset _TMP(PW_15)
set _TMP(PW_16) [pw::BoundaryCondition getByName "bottomwall"]
unset _TMP(PW_16)
set _TMP(PW_17) [pw::BoundaryCondition getByName "top"]
unset _TMP(PW_17)
set _TMP(PW_18) [pw::BoundaryCondition getByName "inlet"]
unset _TMP(PW_18)
set _TMP(PW_19) [pw::BoundaryCondition getByName "outlet"]
unset _TMP(PW_19)
set _TMP(PW_20) [pw::BoundaryCondition getByName "front"]
unset _TMP(PW_20)
set _TMP(PW_21) [pw::BoundaryCondition getByName "back"]
unset _TMP(PW_21)



###############################################################
#-- OUTPUT AND EXIT
#--
###############################################################
pw::Application save [file join $cwd $basename.pw]
set _TMP(mode_2) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1)]]]
  $_TMP(mode_2) initialize -strict -type CAE [file join $cwd $basename.exo]
  $_TMP(mode_2) verify
  $_TMP(mode_2) write
$_TMP(mode_2) end
unset _TMP(mode_2)
