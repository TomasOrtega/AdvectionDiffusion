import Mathlib

set_option linter.unusedVariables false
set_option linter.unusedTactic false
set_option maxHeartbeats 1000000

namespace linear_advection_3d

structure Coordinates where
  x : Real
  y : Real
  z : Real

structure State where
  f : Real

structure Parameters where
  ax : Real
  ay : Real
  az : Real

structure Flux where
  flux_f : Real

structure Fluxes3D where
  xFlux : Flux
  yFlux : Flux
  zFlux : Flux

noncomputable def xFluxExprs (C : Coordinates) (P : Parameters) (U : State) : Flux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    flux_f := (ax * f)
  }

noncomputable def yFluxExprs (C : Coordinates) (P : Parameters) (U : State) : Flux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    flux_f := (ay * f)
  }

noncomputable def zFluxExprs (C : Coordinates) (P : Parameters) (U : State) : Flux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    flux_f := (az * f)
  }

noncomputable def fluxes (C : Coordinates) (P : Parameters) (U : State) : Fluxes3D :=
  {
    xFlux := xFluxExprs C P U
    yFlux := yFluxExprs C P U
    zFlux := zFluxExprs C P U
  }

structure WaveSpeed where
  mu1 : Real

structure WaveSpeedEstimate3D where
  xWaveSpeed : WaveSpeed
  yWaveSpeed : WaveSpeed
  zWaveSpeed : WaveSpeed

noncomputable def xWaveSpeedExprs (C : Coordinates) (P : Parameters) (U : State) : WaveSpeed :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    mu1 := ax
  }

noncomputable def yWaveSpeedExprs (C : Coordinates) (P : Parameters) (U : State) : WaveSpeed :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    mu1 := ay
  }

noncomputable def zWaveSpeedExprs (C : Coordinates) (P : Parameters) (U : State) : WaveSpeed :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    mu1 := az
  }

noncomputable def wavespeeds (C : Coordinates) (P : Parameters) (U : State) : WaveSpeedEstimate3D :=
  {
    xWaveSpeed := xWaveSpeedExprs C P U
    yWaveSpeed := yWaveSpeedExprs C P U
    zWaveSpeed := zWaveSpeedExprs C P U
  }

structure StateGradient where
  d_f : Real

structure FluxJacobian where
  flux_f : StateGradient

structure FluxJacobians3D where
  xFluxJacobian : FluxJacobian
  yFluxJacobian : FluxJacobian
  zFluxJacobian : FluxJacobian

noncomputable def xFluxJacobianExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobian :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    flux_f :=
      {
        d_f := ax
      }
  }

noncomputable def yFluxJacobianExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobian :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    flux_f :=
      {
        d_f := ay
      }
  }

noncomputable def zFluxJacobianExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobian :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    flux_f :=
      {
        d_f := az
      }
  }

noncomputable def fluxJacobians (C : Coordinates) (P : Parameters) (U : State) : FluxJacobians3D :=
  {
    xFluxJacobian := xFluxJacobianExprs C P U
    yFluxJacobian := yFluxJacobianExprs C P U
    zFluxJacobian := zFluxJacobianExprs C P U
  }

structure FluxJacobianEigen where
  lambda1 : Real

structure FluxJacobianEigendecomposition3D where
  xFluxJacobianEigen : FluxJacobianEigen
  yFluxJacobianEigen : FluxJacobianEigen
  zFluxJacobianEigen : FluxJacobianEigen

noncomputable def xFluxJacobianEigenExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigen :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    lambda1 := ax
  }

noncomputable def yFluxJacobianEigenExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigen :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    lambda1 := ay
  }

noncomputable def zFluxJacobianEigenExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigen :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  {
    lambda1 := az
  }

noncomputable def fluxJacobianEigenvalues (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigendecomposition3D :=
  {
    xFluxJacobianEigen := xFluxJacobianEigenExprs C P U
    yFluxJacobianEigen := yFluxJacobianEigenExprs C P U
    zFluxJacobianEigen := zFluxJacobianEigenExprs C P U
  }

theorem xHyperbolicity (C : Coordinates) (P : Parameters) (U : State) :
    (∃ r1 : Real, r1 = (xFluxJacobianEigenExprs C P U).lambda1) := by
  . refine ⟨(xFluxJacobianEigenExprs C P U).lambda1, rfl⟩

theorem yHyperbolicity (C : Coordinates) (P : Parameters) (U : State) :
    (∃ r1 : Real, r1 = (yFluxJacobianEigenExprs C P U).lambda1) := by
  . refine ⟨(yFluxJacobianEigenExprs C P U).lambda1, rfl⟩

theorem zHyperbolicity (C : Coordinates) (P : Parameters) (U : State) :
    (∃ r1 : Real, r1 = (zFluxJacobianEigenExprs C P U).lambda1) := by
  . refine ⟨(zFluxJacobianEigenExprs C P U).lambda1, rfl⟩

theorem xWaveStability (C : Coordinates) (P : Parameters) (U : State) :
    ((abs (xWaveSpeedExprs C P U).mu1) ≥ (abs (xFluxJacobianEigenExprs C P U).lambda1)) := by
  . simp [xWaveSpeedExprs, xFluxJacobianEigenExprs]

theorem yWaveStability (C : Coordinates) (P : Parameters) (U : State) :
    ((abs (yWaveSpeedExprs C P U).mu1) ≥ (abs (yFluxJacobianEigenExprs C P U).lambda1)) := by
  . simp [yWaveSpeedExprs, yFluxJacobianEigenExprs]

theorem zWaveStability (C : Coordinates) (P : Parameters) (U : State) :
    ((abs (zWaveSpeedExprs C P U).mu1) ≥ (abs (zFluxJacobianEigenExprs C P U).lambda1)) := by
  . simp [zWaveSpeedExprs, zFluxJacobianEigenExprs]

structure SpatialGradient where
  f_x : Real
  f_y : Real
  f_z : Real

structure DiffusiveFlux where
  diffusive_flux_f : Real

structure DiffusiveFluxes3D where
  xDiffusiveFlux : DiffusiveFlux
  yDiffusiveFlux : DiffusiveFlux
  zDiffusiveFlux : DiffusiveFlux

noncomputable def xDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  let f_z := DU.f_z
  {
    diffusive_flux_f := 0.0
  }

noncomputable def yDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  let f_z := DU.f_z
  {
    diffusive_flux_f := 0.0
  }

noncomputable def zDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  let f_z := DU.f_z
  {
    diffusive_flux_f := 0.0
  }

noncomputable def diffusiveFluxes (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFluxes3D :=
  {
    xDiffusiveFlux := xDiffusiveFluxExprs C P U DU
    yDiffusiveFlux := yDiffusiveFluxExprs C P U DU
    zDiffusiveFlux := zDiffusiveFluxExprs C P U DU
  }

theorem xDiffusiveFluxConsistency (C : Coordinates) (U : State) (P : Parameters) :
    let DU : SpatialGradient :=
      {
        f_x := 0
        f_y := 0
        f_z := 0
      }
    ((xDiffusiveFluxExprs C P U DU).diffusive_flux_f) = 0 := by
  . simp only [xDiffusiveFluxExprs]
    norm_num

theorem yDiffusiveFluxConsistency (C : Coordinates) (U : State) (P : Parameters) :
    let DU : SpatialGradient :=
      {
        f_x := 0
        f_y := 0
        f_z := 0
      }
    ((yDiffusiveFluxExprs C P U DU).diffusive_flux_f) = 0 := by
  . simp only [yDiffusiveFluxExprs]
    norm_num

theorem zDiffusiveFluxConsistency (C : Coordinates) (U : State) (P : Parameters) :
    let DU : SpatialGradient :=
      {
        f_x := 0
        f_y := 0
        f_z := 0
      }
    ((zDiffusiveFluxExprs C P U DU).diffusive_flux_f) = 0 := by
  . simp only [zDiffusiveFluxExprs]
    norm_num

end linear_advection_3d

namespace linear_advection_3d_lax_friedrichs

open linear_advection_3d

structure Wave where
  wave_f : Real

structure WaveFamily where
  wave1 : Wave
  wave2 : Wave

structure Waves3D where
  xWaveFamily : WaveFamily
  yWaveFamily : WaveFamily
  zWaveFamily : WaveFamily

noncomputable def xWaveFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let z_L := CL.z
  let z_R := CR.z
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    wave1 :=
      {
        wave_f := ((1 / 2 : Real) * ((f_R - f_L) - (((ax * f_R) - (ax * f_L)) / (abs ax))))
      }
    wave2 :=
      {
        wave_f := ((1 / 2 : Real) * ((f_R - f_L) + (((ax * f_R) - (ax * f_L)) / (abs ax))))
      }
  }

noncomputable def yWaveFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let z_L := CL.z
  let z_R := CR.z
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    wave1 :=
      {
        wave_f := ((1 / 2 : Real) * ((f_R - f_L) - (((ay * f_R) - (ay * f_L)) / (abs ay))))
      }
    wave2 :=
      {
        wave_f := ((1 / 2 : Real) * ((f_R - f_L) + (((ay * f_R) - (ay * f_L)) / (abs ay))))
      }
  }

noncomputable def zWaveFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let z_L := CL.z
  let z_R := CR.z
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    wave1 :=
      {
        wave_f := ((1 / 2 : Real) * ((f_R - f_L) - (((az * f_R) - (az * f_L)) / (abs az))))
      }
    wave2 :=
      {
        wave_f := ((1 / 2 : Real) * ((f_R - f_L) + (((az * f_R) - (az * f_L)) / (abs az))))
      }
  }

noncomputable def waves (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Waves3D :=
  {
    xWaveFamily := xWaveFamilyExprs CL CR UL UR P
    yWaveFamily := yWaveFamilyExprs CL CR UL UR P
    zWaveFamily := zWaveFamilyExprs CL CR UL UR P
  }

theorem xWaveConsistency (C : Coordinates) (U : State) (P : Parameters) :
    ((xWaveFamilyExprs C C U U P).wave1.wave_f) = 0.0 ∧
    ((xWaveFamilyExprs C C U U P).wave2.wave_f) = 0.0 := by
  constructor
  . simp [xWaveFamilyExprs]
    ring_nf
    field_simp
  . simp [xWaveFamilyExprs]
    ring_nf
    field_simp

theorem yWaveConsistency (C : Coordinates) (U : State) (P : Parameters) :
    ((yWaveFamilyExprs C C U U P).wave1.wave_f) = 0.0 ∧
    ((yWaveFamilyExprs C C U U P).wave2.wave_f) = 0.0 := by
  constructor
  . simp [yWaveFamilyExprs]
    ring_nf
    field_simp
  . simp [yWaveFamilyExprs]
    ring_nf
    field_simp

theorem zWaveConsistency (C : Coordinates) (U : State) (P : Parameters) :
    ((zWaveFamilyExprs C C U U P).wave1.wave_f) = 0.0 ∧
    ((zWaveFamilyExprs C C U U P).wave2.wave_f) = 0.0 := by
  constructor
  . simp [zWaveFamilyExprs]
    ring_nf
    field_simp
  . simp [zWaveFamilyExprs]
    ring_nf
    field_simp

structure StateJump where
  jump_f : Real

structure StateJumps3D where
  xStateJump : StateJump
  yStateJump : StateJump
  zStateJump : StateJump

noncomputable def xStateJumpExprs (UL : State) (UR : State) : StateJump :=
  let f_L := UL.f
  let f_R := UR.f
  {
    jump_f := f_R - f_L
  }

noncomputable def yStateJumpExprs (UL : State) (UR : State) : StateJump :=
  let f_L := UL.f
  let f_R := UR.f
  {
    jump_f := f_R - f_L
  }

noncomputable def zStateJumpExprs (UL : State) (UR : State) : StateJump :=
  let f_L := UL.f
  let f_R := UR.f
  {
    jump_f := f_R - f_L
  }

noncomputable def stateJumps (UL : State) (UR : State) : StateJumps3D :=
  {
    xStateJump := xStateJumpExprs UL UR
    yStateJump := yStateJumpExprs UL UR
    zStateJump := zStateJumpExprs UL UR
  }

structure WaveSum where
  wave_sum_f : Real

structure WaveSums3D where
  xWaveSum : WaveSum
  yWaveSum : WaveSum
  zWaveSum : WaveSum

noncomputable def xWaveSumExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveSum :=
  {
    wave_sum_f := (xWaveFamilyExprs CL CR UL UR P).wave1.wave_f + (xWaveFamilyExprs CL CR UL UR P).wave2.wave_f
  }

noncomputable def yWaveSumExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveSum :=
  {
    wave_sum_f := (yWaveFamilyExprs CL CR UL UR P).wave1.wave_f + (yWaveFamilyExprs CL CR UL UR P).wave2.wave_f
  }

noncomputable def zWaveSumExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveSum :=
  {
    wave_sum_f := (zWaveFamilyExprs CL CR UL UR P).wave1.wave_f + (zWaveFamilyExprs CL CR UL UR P).wave2.wave_f
  }

noncomputable def waveSums (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveSums3D :=
  {
    xWaveSum := xWaveSumExprs CL CR UL UR P
    yWaveSum := yWaveSumExprs CL CR UL UR P
    zWaveSum := zWaveSumExprs CL CR UL UR P
  }

theorem xWaveJumpCondition (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) :
    ((xStateJumpExprs UL UR).jump_f = (xWaveSumExprs CL CR UL UR P).wave_sum_f) := by
  . simp [xStateJumpExprs, xWaveSumExprs, xWaveFamilyExprs]
    ring_nf

theorem yWaveJumpCondition (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) :
    ((yStateJumpExprs UL UR).jump_f = (yWaveSumExprs CL CR UL UR P).wave_sum_f) := by
  . simp [yStateJumpExprs, yWaveSumExprs, yWaveFamilyExprs]
    ring_nf

theorem zWaveJumpCondition (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) :
    ((zStateJumpExprs UL UR).jump_f = (zWaveSumExprs CL CR UL UR P).wave_sum_f) := by
  . simp [zStateJumpExprs, zWaveSumExprs, zWaveFamilyExprs]
    ring_nf

structure SpeedFamily where
  speed1 : Real
  speed2 : Real

structure Speeds3D where
  xSpeedFamily : SpeedFamily
  ySpeedFamily : SpeedFamily
  zSpeedFamily : SpeedFamily

noncomputable def xSpeedFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : SpeedFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let z_L := CL.z
  let z_R := CR.z
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    speed1 := -((abs ax))
    speed2 := (abs ax)
  }

noncomputable def ySpeedFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : SpeedFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let z_L := CL.z
  let z_R := CR.z
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    speed1 := -((abs ay))
    speed2 := (abs ay)
  }

noncomputable def zSpeedFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : SpeedFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let z_L := CL.z
  let z_R := CR.z
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    speed1 := -((abs az))
    speed2 := (abs az)
  }

noncomputable def speeds (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Speeds3D :=
  {
    xSpeedFamily := xSpeedFamilyExprs CL CR UL UR P
    ySpeedFamily := ySpeedFamilyExprs CL CR UL UR P
    zSpeedFamily := zSpeedFamilyExprs CL CR UL UR P
  }

structure Fluctuation where
  fluctuation_f : Real

structure Fluctuations3D where
  xFluctuation : Fluctuation
  yFluctuation : Fluctuation
  zFluctuation : Fluctuation

noncomputable def xLeftFluctuationExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuation :=
  {
    fluctuation_f := 
      ((xWaveFamilyExprs CL CR UL UR P).wave1.wave_f * (min (xSpeedFamilyExprs CL CR UL UR P).speed1 0.0)) + 
      ((xWaveFamilyExprs CL CR UL UR P).wave2.wave_f * (min (xSpeedFamilyExprs CL CR UL UR P).speed2 0.0))
  }

noncomputable def yLeftFluctuationExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuation :=
  {
    fluctuation_f := 
      ((yWaveFamilyExprs CL CR UL UR P).wave1.wave_f * (min (ySpeedFamilyExprs CL CR UL UR P).speed1 0.0)) + 
      ((yWaveFamilyExprs CL CR UL UR P).wave2.wave_f * (min (ySpeedFamilyExprs CL CR UL UR P).speed2 0.0))
  }

noncomputable def zLeftFluctuationExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuation :=
  {
    fluctuation_f := 
      ((zWaveFamilyExprs CL CR UL UR P).wave1.wave_f * (min (zSpeedFamilyExprs CL CR UL UR P).speed1 0.0)) + 
      ((zWaveFamilyExprs CL CR UL UR P).wave2.wave_f * (min (zSpeedFamilyExprs CL CR UL UR P).speed2 0.0))
  }

noncomputable def xRightFluctuationExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuation :=
  {
    fluctuation_f := 
      ((xWaveFamilyExprs CL CR UL UR P).wave1.wave_f * (max (xSpeedFamilyExprs CL CR UL UR P).speed1 0.0)) + 
      ((xWaveFamilyExprs CL CR UL UR P).wave2.wave_f * (max (xSpeedFamilyExprs CL CR UL UR P).speed2 0.0))
  }

noncomputable def yRightFluctuationExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuation :=
  {
    fluctuation_f := 
      ((yWaveFamilyExprs CL CR UL UR P).wave1.wave_f * (max (ySpeedFamilyExprs CL CR UL UR P).speed1 0.0)) + 
      ((yWaveFamilyExprs CL CR UL UR P).wave2.wave_f * (max (ySpeedFamilyExprs CL CR UL UR P).speed2 0.0))
  }

noncomputable def zRightFluctuationExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuation :=
  {
    fluctuation_f := 
      ((zWaveFamilyExprs CL CR UL UR P).wave1.wave_f * (max (zSpeedFamilyExprs CL CR UL UR P).speed1 0.0)) + 
      ((zWaveFamilyExprs CL CR UL UR P).wave2.wave_f * (max (zSpeedFamilyExprs CL CR UL UR P).speed2 0.0))
  }

noncomputable def leftFluctuations (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuations3D :=
  {
    xFluctuation := xLeftFluctuationExprs CL CR UL UR P
    yFluctuation := yLeftFluctuationExprs CL CR UL UR P
    zFluctuation := zLeftFluctuationExprs CL CR UL UR P
  }

noncomputable def rightFluctuations (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuations3D :=
  {
    xFluctuation := xRightFluctuationExprs CL CR UL UR P
    yFluctuation := yRightFluctuationExprs CL CR UL UR P
    zFluctuation := zRightFluctuationExprs CL CR UL UR P
  }

theorem xLeftFluctuationsConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((xLeftFluctuationExprs C C U U P).fluctuation_f) = 0.0 := by
  . simp [xLeftFluctuationExprs, xWaveFamilyExprs, xSpeedFamilyExprs]
    ring_nf
    field_simp

theorem yLeftFluctuationsConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((yLeftFluctuationExprs C C U U P).fluctuation_f) = 0.0 := by
  . simp [yLeftFluctuationExprs, yWaveFamilyExprs, ySpeedFamilyExprs]
    ring_nf
    field_simp

theorem zLeftFluctuationsConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((zLeftFluctuationExprs C C U U P).fluctuation_f) = 0.0 := by
  . simp [zLeftFluctuationExprs, zWaveFamilyExprs, zSpeedFamilyExprs]
    ring_nf
    field_simp

theorem xRightFluctuationsConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((xRightFluctuationExprs C C U U P).fluctuation_f) = 0.0 := by
  . simp [xRightFluctuationExprs, xWaveFamilyExprs, xSpeedFamilyExprs]
    ring_nf
    field_simp

theorem yRightFluctuationsConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((yRightFluctuationExprs C C U U P).fluctuation_f) = 0.0 := by
  . simp [yRightFluctuationExprs, yWaveFamilyExprs, ySpeedFamilyExprs]
    ring_nf
    field_simp

theorem zRightFluctuationsConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((zRightFluctuationExprs C C U U P).fluctuation_f) = 0.0 := by
  . simp [zRightFluctuationExprs, zWaveFamilyExprs, zSpeedFamilyExprs]
    ring_nf
    field_simp

structure FluxJump where
  flux_jump_f : Real

structure FluxJumps3D where
  xFluxJump : FluxJump
  yFluxJump : FluxJump
  zFluxJump : FluxJump

noncomputable def xFluxJumpExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluxJump :=
  {
    flux_jump_f := (xFluxExprs CR P UR).flux_f - (xFluxExprs CL P UL).flux_f
  }

noncomputable def yFluxJumpExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluxJump :=
  {
    flux_jump_f := (yFluxExprs CR P UR).flux_f - (yFluxExprs CL P UL).flux_f
  }

noncomputable def zFluxJumpExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluxJump :=
  {
    flux_jump_f := (zFluxExprs CR P UR).flux_f - (zFluxExprs CL P UL).flux_f
  }

noncomputable def fluxJumps (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters): FluxJumps3D :=
  {
    xFluxJump := xFluxJumpExprs CL CR UL UR P
    yFluxJump := yFluxJumpExprs CL CR UL UR P
    zFluxJump := zFluxJumpExprs CL CR UL UR P
  }

structure FluctuationSum where
  fluctuation_sum_f : Real

structure FluctuationSums3D where
  xFluctuationSum : FluctuationSum
  yFluctuationSum : FluctuationSum
  zFluctuationSum : FluctuationSum

noncomputable def xFluctuationSumExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluctuationSum :=
  {
    fluctuation_sum_f := (xLeftFluctuationExprs CL CR UL UR P).fluctuation_f + (xRightFluctuationExprs CL CR UL UR P).fluctuation_f
  }

noncomputable def yFluctuationSumExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluctuationSum :=
  {
    fluctuation_sum_f := (yLeftFluctuationExprs CL CR UL UR P).fluctuation_f + (yRightFluctuationExprs CL CR UL UR P).fluctuation_f
  }

noncomputable def zFluctuationSumExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluctuationSum :=
  {
    fluctuation_sum_f := (zLeftFluctuationExprs CL CR UL UR P).fluctuation_f + (zRightFluctuationExprs CL CR UL UR P).fluctuation_f
  }

noncomputable def fluctuationSums (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluctuationSums3D :=
  {
    xFluctuationSum := xFluctuationSumExprs CL CR UL UR P
    yFluctuationSum := yFluctuationSumExprs CL CR UL UR P
    zFluctuationSum := zFluctuationSumExprs CL CR UL UR P
  }

theorem xFluxConservative (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters)
    (Speed1Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let z_L := CL.z
      let z_R := CR.z
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      let az := P.az
      (-((abs ax)) ≠ 0)))
    (Speed2Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let z_L := CL.z
      let z_R := CR.z
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      let az := P.az
      ((abs ax) ≠ 0))) : 
    ((xFluxJumpExprs CL CR UL UR P).flux_jump_f = (xFluctuationSumExprs CL CR UL UR P).fluctuation_sum_f) := by
  . simp [xFluxJumpExprs, xFluctuationSumExprs]
    simp [xFluxExprs, xLeftFluctuationExprs, xRightFluctuationExprs]
    ring_nf
    simp [xWaveFamilyExprs, xSpeedFamilyExprs]
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]
    ring_nf
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]

theorem yFluxConservative (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters)
    (Speed1Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let z_L := CL.z
      let z_R := CR.z
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      let az := P.az
      (-((abs ay)) ≠ 0)))
    (Speed2Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let z_L := CL.z
      let z_R := CR.z
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      let az := P.az
      ((abs ay) ≠ 0))) : 
    ((yFluxJumpExprs CL CR UL UR P).flux_jump_f = (yFluctuationSumExprs CL CR UL UR P).fluctuation_sum_f) := by
  . simp [yFluxJumpExprs, yFluctuationSumExprs]
    simp [yFluxExprs, yLeftFluctuationExprs, yRightFluctuationExprs]
    ring_nf
    simp [yWaveFamilyExprs, ySpeedFamilyExprs]
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]
    ring_nf
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]

theorem zFluxConservative (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters)
    (Speed1Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let z_L := CL.z
      let z_R := CR.z
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      let az := P.az
      (-((abs az)) ≠ 0)))
    (Speed2Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let z_L := CL.z
      let z_R := CR.z
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      let az := P.az
      ((abs az) ≠ 0))) : 
    ((zFluxJumpExprs CL CR UL UR P).flux_jump_f = (zFluctuationSumExprs CL CR UL UR P).fluctuation_sum_f) := by
  . simp [zFluxJumpExprs, zFluctuationSumExprs]
    simp [zFluxExprs, zLeftFluctuationExprs, zRightFluctuationExprs]
    ring_nf
    simp [zWaveFamilyExprs, zSpeedFamilyExprs]
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]
    ring_nf
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]

end linear_advection_3d_lax_friedrichs

namespace linear_advection_3d_minmod

open linear_advection_3d

structure Reconstruction where
  reconstruction_f : Real

structure Reconstruction3D where
  xReconstruction : Reconstruction
  yReconstruction : Reconstruction
  zReconstruction : Reconstruction

noncomputable def xLeftReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let z_L := CL.z
  let z := C.z
  let z_R := CR.z
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    reconstruction_f := (f - ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def yLeftReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let z_L := CL.z
  let z := C.z
  let z_R := CR.z
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    reconstruction_f := (f - ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def zLeftReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let z_L := CL.z
  let z := C.z
  let z_R := CR.z
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    reconstruction_f := (f - ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def xRightReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let z_L := CL.z
  let z := C.z
  let z_R := CR.z
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    reconstruction_f := (f + ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def yRightReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let z_L := CL.z
  let z := C.z
  let z_R := CR.z
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    reconstruction_f := (f + ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def zRightReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let z_L := CL.z
  let z := C.z
  let z_R := CR.z
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  {
    reconstruction_f := (f + ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def leftReconstruction (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction3D :=
  {
    xReconstruction := xLeftReconstructionExprs CL C CR UL U UR P
    yReconstruction := yLeftReconstructionExprs CL C CR UL U UR P
    zReconstruction := zLeftReconstructionExprs CL C CR UL U UR P
  }

noncomputable def rightReconstruction (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction3D :=
  {
    xReconstruction := xRightReconstructionExprs CL C CR UL U UR P
    yReconstruction := yRightReconstructionExprs CL C CR UL U UR P
    zReconstruction := zRightReconstructionExprs CL C CR UL U UR P
  }

theorem xLeftReconstructionConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((xLeftReconstructionExprs C C C U U U P).reconstruction_f) = U.f := by
  . simp [xLeftReconstructionExprs]
    ring_nf
    simp

theorem yLeftReconstructionConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((yLeftReconstructionExprs C C C U U U P).reconstruction_f) = U.f := by
  . simp [yLeftReconstructionExprs]
    ring_nf
    simp

theorem zLeftReconstructionConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((zLeftReconstructionExprs C C C U U U P).reconstruction_f) = U.f := by
  . simp [zLeftReconstructionExprs]
    ring_nf
    simp

theorem xRightReconstructionConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((xRightReconstructionExprs C C C U U U P).reconstruction_f) = U.f := by
  . simp [xRightReconstructionExprs]
    ring_nf
    simp

theorem yRightReconstructionConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((yRightReconstructionExprs C C C U U U P).reconstruction_f) = U.f := by
  . simp [yRightReconstructionExprs]
    ring_nf
    simp

theorem zRightReconstructionConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((zRightReconstructionExprs C C C U U U P).reconstruction_f) = U.f := by
  . simp [zRightReconstructionExprs]
    ring_nf
    simp

theorem xLeftReconstructionLinearityPreservation (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (U : State) (P : Parameters) (D : Real) :
    let UL : State :=
      {
        f := U.f - D
      }
    let UR : State :=
      {
        f := U.f + D
      }
    ((xLeftReconstructionExprs CL C CR UL U UR P).reconstruction_f) = U.f - ((1 / 2 : Real) * D) := by
  . simp [xLeftReconstructionExprs]
    ring_nf
    simp

theorem yLeftReconstructionLinearityPreservation (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (U : State) (P : Parameters) (D : Real) :
    let UL : State :=
      {
        f := U.f - D
      }
    let UR : State :=
      {
        f := U.f + D
      }
    ((yLeftReconstructionExprs CL C CR UL U UR P).reconstruction_f) = U.f - ((1 / 2 : Real) * D) := by
  . simp [yLeftReconstructionExprs]
    ring_nf
    simp

theorem zLeftReconstructionLinearityPreservation (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (U : State) (P : Parameters) (D : Real) :
    let UL : State :=
      {
        f := U.f - D
      }
    let UR : State :=
      {
        f := U.f + D
      }
    ((zLeftReconstructionExprs CL C CR UL U UR P).reconstruction_f) = U.f - ((1 / 2 : Real) * D) := by
  . simp [zLeftReconstructionExprs]
    ring_nf
    simp

theorem xRightReconstructionLinearityPreservation (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (U : State) (P : Parameters) (D : Real) :
    let UL : State :=
      {
        f := U.f - D
      }
    let UR : State :=
      {
        f := U.f + D
      }
    ((xRightReconstructionExprs CL C CR UL U UR P).reconstruction_f) = U.f + ((1 / 2 : Real) * D) := by
  . simp [xRightReconstructionExprs]
    ring_nf
    simp

theorem yRightReconstructionLinearityPreservation (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (U : State) (P : Parameters) (D : Real) :
    let UL : State :=
      {
        f := U.f - D
      }
    let UR : State :=
      {
        f := U.f + D
      }
    ((yRightReconstructionExprs CL C CR UL U UR P).reconstruction_f) = U.f + ((1 / 2 : Real) * D) := by
  . simp [yRightReconstructionExprs]
    ring_nf
    simp

theorem zRightReconstructionLinearityPreservation (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (U : State) (P : Parameters) (D : Real) :
    let UL : State :=
      {
        f := U.f - D
      }
    let UR : State :=
      {
        f := U.f + D
      }
    ((zRightReconstructionExprs CL C CR UL U UR P).reconstruction_f) = U.f + ((1 / 2 : Real) * D) := by
  . simp [zRightReconstructionExprs]
    ring_nf
    simp

theorem xReconstructionSymmetric (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) :
    ((xLeftReconstructionExprs CL C CR UL U UR P).reconstruction_f) = ((xRightReconstructionExprs CR C CL UR U UL P).reconstruction_f) := by
  . simp [xLeftReconstructionExprs, xRightReconstructionExprs]
    ring_nf
    field_simp
    simp only [min_def, max_def]
    split_ifs <;> linarith

theorem yReconstructionSymmetric (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) :
    ((yLeftReconstructionExprs CL C CR UL U UR P).reconstruction_f) = ((yRightReconstructionExprs CR C CL UR U UL P).reconstruction_f) := by
  . simp [yLeftReconstructionExprs, yRightReconstructionExprs]
    ring_nf
    field_simp
    simp only [min_def, max_def]
    split_ifs <;> linarith

theorem zReconstructionSymmetric (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) :
    ((zLeftReconstructionExprs CL C CR UL U UR P).reconstruction_f) = ((zRightReconstructionExprs CR C CL UR U UL P).reconstruction_f) := by
  . simp [zLeftReconstructionExprs, zRightReconstructionExprs]
    ring_nf
    field_simp
    simp only [min_def, max_def]
    split_ifs <;> linarith

end linear_advection_3d_minmod
                