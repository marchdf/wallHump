# Pointwise V18.0R4 Journal file - Tue Nov 14 19:00:09 2017

package require PWI_Glyph 2.18.0

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

pw::Application reset -keep Clipboard
set _TMP(mode_2) [pw::Application begin ProjectLoader]
  $_TMP(mode_2) initialize {/Users/mhenryde/wind/wallHump/threed/103x28/hump2newtop_noplenumZ103x28_3D.pw}
  $_TMP(mode_2) setAppendMode false
  $_TMP(mode_2) setRepairMode Defer
  $_TMP(mode_2) load
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application resetUndoLevels
set _TMP(mode_3) [pw::Application begin Create]
  set _DM(1) [pw::GridEntity getByName "dom-1"]
  set _TMP(PW_1) [pw::FaceStructured createFromDomains [list $_DM(1)]]
  set _TMP(face_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _TMP(extStrBlock_1) [pw::BlockStructured create]
  $_TMP(extStrBlock_1) addFace $_TMP(face_1)
$_TMP(mode_3) end
unset _TMP(mode_3)
set _TMP(mode_4) [pw::Application begin ExtrusionSolver [list $_TMP(extStrBlock_1)]]
  $_TMP(mode_4) setKeepFailingStep true
  $_TMP(extStrBlock_1) setExtrusionSolverAttribute Mode Translate
  $_TMP(extStrBlock_1) setExtrusionSolverAttribute TranslateDirection {1 0 0}
  set _BL(1) [pw::GridEntity getByName "blk-1"]
  $_TMP(extStrBlock_1) setExtrusionSolverAttribute TranslateDirection {0 0 1}
  $_TMP(extStrBlock_1) setExtrusionSolverAttribute TranslateDistance 0.25
  $_TMP(mode_4) run 25
$_TMP(mode_4) end
unset _TMP(mode_4)
unset _TMP(extStrBlock_1)
unset _TMP(face_1)
pw::Application markUndoLevel {Extrude, Translate}

set _DM(2) [pw::GridEntity getByName "dom-2"]
set _DM(3) [pw::GridEntity getByName "dom-3"]
set _DM(4) [pw::GridEntity getByName "dom-4"]
set _DM(5) [pw::GridEntity getByName "dom-5"]
set _DM(6) [pw::GridEntity getByName "dom-6"]
set _TMP(PW_2) [pw::BoundaryCondition getByName {Unspecified}]
set _TMP(PW_3) [pw::BoundaryCondition getByName {bottomwall}]
set _TMP(PW_4) [pw::BoundaryCondition getByName {top}]
set _TMP(PW_5) [pw::BoundaryCondition getByName {inlet}]
set _TMP(PW_6) [pw::BoundaryCondition getByName {outlet}]
set _TMP(PW_7) [pw::BoundaryCondition getByName {front}]
set _TMP(PW_8) [pw::BoundaryCondition getByName {back}]
$_TMP(PW_5) apply [list [list $_BL(1) $_DM(5)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_4) apply [list [list $_BL(1) $_DM(4)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_2) apply [list [list $_BL(1) $_DM(1)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_7) apply [list [list $_BL(1) $_DM(6)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_8) apply [list [list $_BL(1) $_DM(1)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_6) apply [list [list $_BL(1) $_DM(3)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_3) apply [list [list $_BL(1) $_DM(2)]]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_2)
unset _TMP(PW_3)
unset _TMP(PW_4)
unset _TMP(PW_5)
unset _TMP(PW_6)
unset _TMP(PW_7)
unset _TMP(PW_8)
set _TMP(mode_5) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1)]]]
  $_TMP(mode_5) initialize -strict -type CAE {/Users/mhenryde/wind/wallHump/threed/103x28/hump2newtop_noplenumZ103x28x25_3D.exo}
  $_TMP(mode_5) verify
  $_TMP(mode_5) write
$_TMP(mode_5) end
unset _TMP(mode_5)
