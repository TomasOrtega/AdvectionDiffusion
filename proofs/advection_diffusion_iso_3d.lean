import Mathlib

set_option linter.unusedVariables false
set_option linter.unusedTactic false
set_option maxHeartbeats 1000000

namespace advection_diffusion_iso_3d

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
  D : Real

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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  let f_z := DU.f_z
  {
    diffusive_flux_f := (D * f_x)
  }

noncomputable def yDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  let f_z := DU.f_z
  {
    diffusive_flux_f := (D * f_y)
  }

noncomputable def zDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let z := C.z
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  let f_z := DU.f_z
  {
    diffusive_flux_f := (D * f_z)
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

end advection_diffusion_iso_3d

namespace advection_diffusion_iso_3d_lax_friedrichs

open advection_diffusion_iso_3d

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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
      let D := P.D
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
      let D := P.D
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
      let D := P.D
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
      let D := P.D
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
      let D := P.D
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
      let D := P.D
      ((abs az) ≠ 0))) : 
    ((zFluxJumpExprs CL CR UL UR P).flux_jump_f = (zFluctuationSumExprs CL CR UL UR P).fluctuation_sum_f) := by
  . simp [zFluxJumpExprs, zFluctuationSumExprs]
    simp [zFluxExprs, zLeftFluctuationExprs, zRightFluctuationExprs]
    ring_nf
    simp [zWaveFamilyExprs, zSpeedFamilyExprs]
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]
    ring_nf
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]

end advection_diffusion_iso_3d_lax_friedrichs

namespace advection_diffusion_iso_3d_minmod

open advection_diffusion_iso_3d

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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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
  let D := P.D
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

end advection_diffusion_iso_3d_minmod

namespace advection_diffusion_iso_3d_centered_diffusion

open advection_diffusion_iso_3d

structure InterfaceState where
  interface_f : Real

structure InterfaceState3D where
  xInterfaceState : InterfaceState
  yInterfaceState : InterfaceState
  zInterfaceState : InterfaceState

noncomputable def xInterfaceStateExprs (C_Lx : Coordinates) (C_Rx : Coordinates)
    (U_Lx : State) (U_Rx : State) (P : Parameters) : InterfaceState :=
  let x_Lx := C_Lx.x
  let x_Rx := C_Rx.x
  let y_Lx := C_Lx.y
  let y_Rx := C_Rx.y
  let z_Lx := C_Lx.z
  let z_Rx := C_Rx.z
  let f_Lx := U_Lx.f
  let f_Rx := U_Rx.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  {
    interface_f := ((1 / 2 : Real) * (f_Lx + f_Rx))
  }

noncomputable def yInterfaceStateExprs (C_Ly : Coordinates) (C_Ry : Coordinates)
    (U_Ly : State) (U_Ry : State) (P : Parameters) : InterfaceState :=
  let x_Ly := C_Ly.x
  let x_Ry := C_Ry.x
  let y_Ly := C_Ly.y
  let y_Ry := C_Ry.y
  let z_Ly := C_Ly.z
  let z_Ry := C_Ry.z
  let f_Ly := U_Ly.f
  let f_Ry := U_Ry.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  {
    interface_f := ((1 / 2 : Real) * (f_Ly + f_Ry))
  }

noncomputable def zInterfaceStateExprs (C_Lz : Coordinates) (C_Rz : Coordinates)
    (U_Lz : State) (U_Rz : State) (P : Parameters) : InterfaceState :=
  let x_Lz := C_Lz.x
  let x_Rz := C_Rz.x
  let y_Lz := C_Lz.y
  let y_Rz := C_Rz.y
  let z_Lz := C_Lz.z
  let z_Rz := C_Rz.z
  let f_Lz := U_Lz.f
  let f_Rz := U_Rz.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  {
    interface_f := ((1 / 2 : Real) * (f_Lz + f_Rz))
  }

noncomputable def interfaceState (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : InterfaceState3D :=
  {
    xInterfaceState := xInterfaceStateExprs CL CR UL UR P
    yInterfaceState := yInterfaceStateExprs CL CR UL UR P
    zInterfaceState := zInterfaceStateExprs CL CR UL UR P
  }

theorem xInterfaceStateConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((xInterfaceStateExprs C C U U P).interface_f) = U.f := by
  . simp [xInterfaceStateExprs]
    ring_nf

theorem yInterfaceStateConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((yInterfaceStateExprs C C U U P).interface_f) = U.f := by
  . simp [yInterfaceStateExprs]
    ring_nf

theorem zInterfaceStateConsistent (C : Coordinates) (U : State) (P : Parameters) :
    ((zInterfaceStateExprs C C U U P).interface_f) = U.f := by
  . simp [zInterfaceStateExprs]
    ring_nf

theorem xInterfaceStateSymmetric (C_Lx : Coordinates) (C_Rx : Coordinates)
    (U_Lx : State) (U_Rx : State) (P : Parameters) :
    ((xInterfaceStateExprs C_Lx C_Rx U_Lx U_Rx P).interface_f) = ((xInterfaceStateExprs C_Rx C_Lx U_Rx U_Lx P).interface_f) := by
  . simp [xInterfaceStateExprs]
    ring_nf

theorem yInterfaceStateSymmetric (C_Ly : Coordinates) (C_Ry : Coordinates)
    (U_Ly : State) (U_Ry : State) (P : Parameters) :
    ((yInterfaceStateExprs C_Ly C_Ry U_Ly U_Ry P).interface_f) = ((yInterfaceStateExprs C_Ry C_Ly U_Ry U_Ly P).interface_f) := by
  . simp [yInterfaceStateExprs]
    ring_nf

theorem zInterfaceStateSymmetric (C_Lz : Coordinates) (C_Rz : Coordinates)
    (U_Lz : State) (U_Rz : State) (P : Parameters) :
    ((zInterfaceStateExprs C_Lz C_Rz U_Lz U_Rz P).interface_f) = ((zInterfaceStateExprs C_Rz C_Lz U_Rz U_Lz P).interface_f) := by
  . simp [zInterfaceStateExprs]
    ring_nf

structure InterfaceGradient where
  interface_f_x : Real
  interface_f_y : Real
  interface_f_z : Real

structure InterfaceGradients3D where
  xInterfaceGradient : InterfaceGradient
  yInterfaceGradient : InterfaceGradient
  zInterfaceGradient : InterfaceGradient

noncomputable def xInterfaceGradientExprs (C_Lx : Coordinates) (C_Rx : Coordinates)
    (U_Lx: State) (U_Rx : State)
    (C_Lx_Ly : Coordinates) (C_Lx_Ry : Coordinates)
    (C_Rx_Ly : Coordinates) (C_Rx_Ry : Coordinates)
    (U_Lx_Ly : State) (U_Lx_Ry : State)
    (U_Rx_Ly : State) (U_Rx_Ry : State)
    (C_Lx_Lz : Coordinates) (C_Lx_Rz : Coordinates)
    (C_Rx_Lz : Coordinates) (C_Rx_Rz : Coordinates)
    (U_Lx_Lz : State) (U_Lx_Rz : State)
    (U_Rx_Lz : State) (U_Rx_Rz : State) (P : Parameters) : InterfaceGradient :=
  let x_Lx := C_Lx.x
  let x_Rx := C_Rx.x
  let y_Lx_Ly := C_Lx_Ly.y
  let y_Lx_Ry := C_Lx_Ry.y
  let y_Rx_Ly := C_Rx_Ly.y
  let y_Rx_Ry := C_Rx_Ry.y
  let z_Lx_Lz := C_Lx_Lz.z
  let z_Lx_Rz := C_Lx_Rz.z
  let z_Rx_Lz := C_Rx_Lz.z
  let z_Rx_Rz := C_Rx_Rz.z
  let f_Lx := U_Lx.f
  let f_Rx := U_Rx.f
  let f_Lx_Ly := U_Lx_Ly.f
  let f_Lx_Ry := U_Lx_Ry.f
  let f_Rx_Ly := U_Rx_Ly.f
  let f_Rx_Ry := U_Rx_Ry.f
  let f_Lx_Lz := U_Lx_Lz.f
  let f_Lx_Rz := U_Lx_Rz.f
  let f_Rx_Lz := U_Rx_Lz.f
  let f_Rx_Rz := U_Rx_Rz.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  {
    interface_f_x := ((f_Rx - f_Lx) / (x_Rx - x_Lx))
    interface_f_y := ((1 / 2 : Real) * (((f_Lx_Ry - f_Lx_Ly) / (y_Lx_Ry - y_Lx_Ly)) + ((f_Rx_Ry - f_Rx_Ly) / (y_Rx_Ry -
      y_Rx_Ly))))
    interface_f_z := ((1 / 2 : Real) * (((f_Lx_Rz - f_Lx_Lz) / (z_Lx_Rz - z_Lx_Lz)) + ((f_Rx_Rz - f_Rx_Lz) / (z_Rx_Rz -
      z_Rx_Lz))))
  }

noncomputable def yInterfaceGradientExprs (C_Ly_Lx : Coordinates) (C_Ly_Rx : Coordinates)
    (C_Ry_Lx : Coordinates) (C_Ry_Rx : Coordinates)
    (U_Ly_Lx : State) (U_Ly_Rx : State)
    (U_Ry_Lx : State) (U_Ry_Rx : State)
    (C_Ly : Coordinates) (C_Ry : Coordinates)
    (U_Ly: State) (U_Ry : State)
    (C_Ly_Lz : Coordinates) (C_Ly_Rz : Coordinates)
    (C_Ry_Lz : Coordinates) (C_Ry_Rz : Coordinates)
    (U_Ly_Lz : State) (U_Ly_Rz : State)
    (U_Ry_Lz : State) (U_Ry_Rz : State) (P : Parameters) : InterfaceGradient :=
  let x_Ly_Lx := C_Ly_Lx.x
  let x_Ly_Rx := C_Ly_Rx.x
  let x_Ry_Lx := C_Ry_Lx.x
  let x_Ry_Rx := C_Ry_Rx.x
  let y_Ly := C_Ly.y
  let y_Ry := C_Ry.y
  let z_Ly_Lz := C_Ly_Lz.z
  let z_Ly_Rz := C_Ly_Rz.z
  let z_Ry_Lz := C_Ry_Lz.z
  let z_Ry_Rz := C_Ry_Rz.z
  let f_Ly_Lx := U_Ly_Lx.f
  let f_Ly_Rx := U_Ly_Rx.f
  let f_Ry_Lx := U_Ry_Lx.f
  let f_Ry_Rx := U_Ry_Rx.f
  let f_Ly := U_Ly.f
  let f_Ry := U_Ry.f
  let f_Ly_Lz := U_Ly_Lz.f
  let f_Ly_Rz := U_Ly_Rz.f
  let f_Ry_Lz := U_Ry_Lz.f
  let f_Ry_Rz := U_Ry_Rz.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  {
    interface_f_x := ((1 / 2 : Real) * (((f_Ly_Rx - f_Ly_Lx) / (x_Ly_Rx - x_Ly_Lx)) + ((f_Ry_Rx - f_Ry_Lx) / (x_Ry_Rx -
      x_Ry_Lx))))
    interface_f_y := ((f_Ry - f_Ly) / (y_Ry - y_Ly))
    interface_f_z := ((1 / 2 : Real) * (((f_Ly_Rz - f_Ly_Lz) / (z_Ly_Rz - z_Ly_Lz)) + ((f_Ry_Rz - f_Ry_Lz) / (z_Ry_Rz -
      z_Ry_Lz))))
  }

noncomputable def zInterfaceGradientExprs (C_Lz_Lx : Coordinates) (C_Lz_Rx : Coordinates)
    (C_Rz_Lx : Coordinates) (C_Rz_Rx : Coordinates)
    (U_Lz_Lx : State) (U_Lz_Rx : State)
    (U_Rz_Lx : State) (U_Rz_Rx : State)
    (C_Lz_Ly : Coordinates) (C_Lz_Ry : Coordinates)
    (C_Rz_Ly : Coordinates) (C_Rz_Ry : Coordinates)
    (U_Lz_Ly : State) (U_Lz_Ry : State)
    (U_Rz_Ly : State) (U_Rz_Ry : State)
    (C_Lz : Coordinates) (C_Rz : Coordinates)
    (U_Lz: State) (U_Rz : State) (P : Parameters) : InterfaceGradient :=
  let x_Lz_Lx := C_Lz_Lx.x
  let x_Lz_Rx := C_Lz_Rx.x
  let x_Rz_Lx := C_Rz_Lx.x
  let x_Rz_Rx := C_Rz_Rx.x
  let y_Lz_Ly := C_Lz_Ly.y
  let y_Lz_Ry := C_Lz_Ry.y
  let y_Rz_Ly := C_Rz_Ly.y
  let y_Rz_Ry := C_Rz_Ry.y
  let z_Lz := C_Lz.z
  let z_Rz := C_Rz.z
  let f_Lz_Lx := U_Lz_Lx.f
  let f_Lz_Rx := U_Lz_Rx.f
  let f_Rz_Lx := U_Rz_Lx.f
  let f_Rz_Rx := U_Rz_Rx.f
  let f_Lz_Ly := U_Lz_Ly.f
  let f_Lz_Ry := U_Lz_Ry.f
  let f_Rz_Ly := U_Rz_Ly.f
  let f_Rz_Ry := U_Rz_Ry.f
  let f_Lz := U_Lz.f
  let f_Rz := U_Rz.f
  let ax := P.ax
  let ay := P.ay
  let az := P.az
  let D := P.D
  {
    interface_f_x := ((1 / 2 : Real) * (((f_Lz_Rx - f_Lz_Lx) / (x_Lz_Rx - x_Lz_Lx)) + ((f_Rz_Rx - f_Rz_Lx) / (x_Rz_Rx -
      x_Rz_Lx))))
    interface_f_y := ((1 / 2 : Real) * (((f_Lz_Ry - f_Lz_Ly) / (y_Lz_Ry - y_Lz_Ly)) + ((f_Rz_Ry - f_Rz_Ly) / (y_Rz_Ry -
      y_Rz_Ly))))
    interface_f_z := ((f_Rz - f_Lz) / (z_Rz - z_Lz))
  }

theorem xInterfaceGradientConsistency (C_Lx : Coordinates) (C_Rx : Coordinates)
    (Ux: State)
    (C_Lx_Ly : Coordinates) (C_Lx_Ry : Coordinates)
    (C_Rx_Ly : Coordinates) (C_Rx_Ry : Coordinates)
    (U_Lx : State) (U_Rx : State)
    (C_Lx_Lz : Coordinates) (C_Lx_Rz : Coordinates)
    (C_Rx_Lz : Coordinates) (C_Rx_Rz : Coordinates)
    (U_Lx : State) (U_Rx : State) (P : Parameters) :
    ((xInterfaceGradientExprs C_Lx C_Rx Ux Ux C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx U_Lx U_Rx U_Rx C_Lx_Lz C_Lx_Rz
      C_Rx_Lz C_Rx_Rz U_Lx U_Lx U_Rx U_Rx P).interface_f_x) = 0.0 ∧
    ((xInterfaceGradientExprs C_Lx C_Rx Ux Ux C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx U_Lx U_Rx U_Rx C_Lx_Lz C_Lx_Rz
      C_Rx_Lz C_Rx_Rz U_Lx U_Lx U_Rx U_Rx P).interface_f_y) = 0.0 ∧
    ((xInterfaceGradientExprs C_Lx C_Rx Ux Ux C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx U_Lx U_Rx U_Rx C_Lx_Lz C_Lx_Rz
      C_Rx_Lz C_Rx_Rz U_Lx U_Lx U_Rx U_Rx P).interface_f_z) = 0.0 := by
  constructor
  . simp [xInterfaceGradientExprs]
    ring_nf
    field_simp
  constructor
  . simp [xInterfaceGradientExprs]
    ring_nf
    field_simp
  . simp [xInterfaceGradientExprs]
    ring_nf
    field_simp

theorem yInterfaceGradientConsistency (C_Ly_Lx : Coordinates) (C_Ly_Rx : Coordinates)
    (C_Ry_Lx : Coordinates) (C_Ry_Rx : Coordinates)
    (U_Ly : State) (U_Ry : State)
    (C_Ly : Coordinates) (C_Ry : Coordinates)
    (Uy: State)
    (C_Ly_Lz : Coordinates) (C_Ly_Rz : Coordinates)
    (C_Ry_Lz : Coordinates) (C_Ry_Rz : Coordinates)
    (U_Ly : State) (U_Ry : State) (P : Parameters) :
    ((yInterfaceGradientExprs C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly U_Ly U_Ry U_Ry C_Ly C_Ry Uy Uy C_Ly_Lz C_Ly_Rz
      C_Ry_Lz C_Ry_Rz U_Ly U_Ly U_Ry U_Ry P).interface_f_x) = 0.0 ∧
    ((yInterfaceGradientExprs C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly U_Ly U_Ry U_Ry C_Ly C_Ry Uy Uy C_Ly_Lz C_Ly_Rz
      C_Ry_Lz C_Ry_Rz U_Ly U_Ly U_Ry U_Ry P).interface_f_y) = 0.0 ∧
    ((yInterfaceGradientExprs C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly U_Ly U_Ry U_Ry C_Ly C_Ry Uy Uy C_Ly_Lz C_Ly_Rz
      C_Ry_Lz C_Ry_Rz U_Ly U_Ly U_Ry U_Ry P).interface_f_z) = 0.0 := by
  constructor
  . simp [yInterfaceGradientExprs]
    ring_nf
    field_simp
  constructor
  . simp [yInterfaceGradientExprs]
    ring_nf
    field_simp
  . simp [yInterfaceGradientExprs]
    ring_nf
    field_simp

theorem zInterfaceGradientConsistency (C_Lz_Lx : Coordinates) (C_Lz_Rx : Coordinates)
    (C_Rz_Lx : Coordinates) (C_Rz_Rx : Coordinates)
    (U_Lz : State) (U_Rz : State)
    (C_Lz_Ly : Coordinates) (C_Lz_Ry : Coordinates)
    (C_Rz_Ly : Coordinates) (C_Rz_Ry : Coordinates)
    (U_Lz : State) (U_Rz : State)
    (C_Lz : Coordinates) (C_Rz : Coordinates)
    (Uz: State) (P : Parameters) :
    ((zInterfaceGradientExprs C_Lz_Lx C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz U_Lz U_Rz U_Rz C_Lz_Ly C_Lz_Ry C_Rz_Ly C_Rz_Ry U_Lz
      U_Lz U_Rz U_Rz C_Lz C_Rz Uz Uz P).interface_f_x) = 0.0 ∧
    ((zInterfaceGradientExprs C_Lz_Lx C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz U_Lz U_Rz U_Rz C_Lz_Ly C_Lz_Ry C_Rz_Ly C_Rz_Ry U_Lz
      U_Lz U_Rz U_Rz C_Lz C_Rz Uz Uz P).interface_f_y) = 0.0 ∧
    ((zInterfaceGradientExprs C_Lz_Lx C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz U_Lz U_Rz U_Rz C_Lz_Ly C_Lz_Ry C_Rz_Ly C_Rz_Ry U_Lz
      U_Lz U_Rz U_Rz C_Lz C_Rz Uz Uz P).interface_f_z) = 0.0 := by
  constructor
  . simp [zInterfaceGradientExprs]
    ring_nf
    field_simp
  constructor
  . simp [zInterfaceGradientExprs]
    ring_nf
    field_simp
  . simp [zInterfaceGradientExprs]
    ring_nf
    field_simp

theorem xInterfaceGradientAntisymmetric (C_Lx : Coordinates) (C_Rx : Coordinates)
    (U_Lx : State) (U_Rx : State)
    (C_Lx_Ly : Coordinates) (C_Lx_Ry : Coordinates)
    (C_Rx_Ly : Coordinates) (C_Rx_Ry : Coordinates)
    (U_Lx_Ly : State) (U_Lx_Ry : State)
    (U_Rx_Ly : State) (U_Rx_Ry : State)
    (C_Lx_Lz : Coordinates) (C_Lx_Rz : Coordinates)
    (C_Rx_Lz : Coordinates) (C_Rx_Rz : Coordinates)
    (U_Lx_Lz : State) (U_Lx_Rz : State)
    (U_Rx_Lz : State) (U_Rx_Rz : State) (P : Parameters) :
    ((xInterfaceGradientExprs C_Lx C_Rx U_Lx U_Rx C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx_Ly U_Lx_Ry U_Rx_Ly U_Rx_Ry
      C_Lx_Lz C_Lx_Rz C_Rx_Lz C_Rx_Rz U_Lx_Lz U_Lx_Rz U_Rx_Lz U_Rx_Rz P).interface_f_x) = -((xInterfaceGradientExprs
      C_Lx C_Rx U_Rx U_Lx C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx_Ry U_Lx_Ly U_Rx_Ry U_Rx_Ly C_Lx_Lz C_Lx_Rz C_Rx_Lz
      C_Rx_Rz U_Lx_Rz U_Lx_Lz U_Rx_Rz U_Rx_Lz P).interface_f_x) ∧
    ((xInterfaceGradientExprs C_Lx C_Rx U_Lx U_Rx C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx_Ly U_Lx_Ry U_Rx_Ly U_Rx_Ry
      C_Lx_Lz C_Lx_Rz C_Rx_Lz C_Rx_Rz U_Lx_Lz U_Lx_Rz U_Rx_Lz U_Rx_Rz P).interface_f_y) = -((xInterfaceGradientExprs
      C_Lx C_Rx U_Rx U_Lx C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx_Ry U_Lx_Ly U_Rx_Ry U_Rx_Ly C_Lx_Lz C_Lx_Rz C_Rx_Lz
      C_Rx_Rz U_Lx_Rz U_Lx_Lz U_Rx_Rz U_Rx_Lz P).interface_f_y) ∧
    ((xInterfaceGradientExprs C_Lx C_Rx U_Lx U_Rx C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx_Ly U_Lx_Ry U_Rx_Ly U_Rx_Ry
      C_Lx_Lz C_Lx_Rz C_Rx_Lz C_Rx_Rz U_Lx_Lz U_Lx_Rz U_Rx_Lz U_Rx_Rz P).interface_f_z) = -((xInterfaceGradientExprs
      C_Lx C_Rx U_Rx U_Lx C_Lx_Ly C_Lx_Ry C_Rx_Ly C_Rx_Ry U_Lx_Ry U_Lx_Ly U_Rx_Ry U_Rx_Ly C_Lx_Lz C_Lx_Rz C_Rx_Lz
      C_Rx_Rz U_Lx_Rz U_Lx_Lz U_Rx_Rz U_Rx_Lz P).interface_f_z) := by
  constructor
  . simp [xInterfaceGradientExprs]
    ring_nf
  constructor
  . simp [xInterfaceGradientExprs]
    ring_nf
  . simp [xInterfaceGradientExprs]
    ring_nf

theorem yInterfaceGradientAntisymmetric (C_Ly_Lx : Coordinates) (C_Ly_Rx : Coordinates)
    (C_Ry_Lx : Coordinates) (C_Ry_Rx : Coordinates)
    (U_Ly_Lx : State) (U_Ly_Rx : State)
    (U_Ry_Lx : State) (U_Ry_Rx : State)
    (C_Ly : Coordinates) (C_Ry : Coordinates)
    (U_Ly : State) (U_Ry : State)
    (C_Ly_Lz : Coordinates) (C_Ly_Rz : Coordinates)
    (C_Ry_Lz : Coordinates) (C_Ry_Rz : Coordinates)
    (U_Ly_Lz : State) (U_Ly_Rz : State)
    (U_Ry_Lz : State) (U_Ry_Rz : State) (P : Parameters) :
    ((yInterfaceGradientExprs C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly_Lx U_Ly_Rx U_Ry_Lx U_Ry_Rx C_Ly C_Ry U_Ly U_Ry
      C_Ly_Lz C_Ly_Rz C_Ry_Lz C_Ry_Rz U_Ly_Lz U_Ly_Rz U_Ry_Lz U_Ry_Rz P).interface_f_x) = -((yInterfaceGradientExprs
      C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly_Rx U_Ly_Lx U_Ry_Rx U_Ry_Lx C_Ly C_Ry U_Ry U_Ly C_Ly_Lz C_Ly_Rz C_Ry_Lz
      C_Ry_Rz U_Ly_Rz U_Ly_Lz U_Ry_Rz U_Ry_Lz P).interface_f_x) ∧
    ((yInterfaceGradientExprs C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly_Lx U_Ly_Rx U_Ry_Lx U_Ry_Rx C_Ly C_Ry U_Ly U_Ry
      C_Ly_Lz C_Ly_Rz C_Ry_Lz C_Ry_Rz U_Ly_Lz U_Ly_Rz U_Ry_Lz U_Ry_Rz P).interface_f_y) = -((yInterfaceGradientExprs
      C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly_Rx U_Ly_Lx U_Ry_Rx U_Ry_Lx C_Ly C_Ry U_Ry U_Ly C_Ly_Lz C_Ly_Rz C_Ry_Lz
      C_Ry_Rz U_Ly_Rz U_Ly_Lz U_Ry_Rz U_Ry_Lz P).interface_f_y) ∧
    ((yInterfaceGradientExprs C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly_Lx U_Ly_Rx U_Ry_Lx U_Ry_Rx C_Ly C_Ry U_Ly U_Ry
      C_Ly_Lz C_Ly_Rz C_Ry_Lz C_Ry_Rz U_Ly_Lz U_Ly_Rz U_Ry_Lz U_Ry_Rz P).interface_f_z) = -((yInterfaceGradientExprs
      C_Ly_Lx C_Ly_Rx C_Ry_Lx C_Ry_Rx U_Ly_Rx U_Ly_Lx U_Ry_Rx U_Ry_Lx C_Ly C_Ry U_Ry U_Ly C_Ly_Lz C_Ly_Rz C_Ry_Lz
      C_Ry_Rz U_Ly_Rz U_Ly_Lz U_Ry_Rz U_Ry_Lz P).interface_f_z) := by
  constructor
  . simp [yInterfaceGradientExprs]
    ring_nf
  constructor
  . simp [yInterfaceGradientExprs]
    ring_nf
  . simp [yInterfaceGradientExprs]
    ring_nf

theorem zInterfaceGradientAntisymmetric (C_Lz_Lx : Coordinates) (C_Lz_Rx : Coordinates)
    (C_Rz_Lx : Coordinates) (C_Rz_Rx : Coordinates)
    (U_Lz_Lx : State) (U_Lz_Rx : State)
    (U_Rz_Lx : State) (U_Rz_Rx : State)
    (C_Lz_Ly : Coordinates) (C_Lz_Ry : Coordinates)
    (C_Rz_Ly : Coordinates) (C_Rz_Ry : Coordinates)
    (U_Lz_Ly : State) (U_Lz_Ry : State)
    (U_Rz_Ly : State) (U_Rz_Ry : State)
    (C_Lz : Coordinates) (C_Rz : Coordinates)
    (U_Lz : State) (U_Rz : State) (P : Parameters) :
    ((zInterfaceGradientExprs C_Lz_Lx C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz_Lx U_Lz_Rx U_Rz_Lx U_Rz_Rx C_Lz_Ly C_Lz_Ry C_Rz_Ly
      C_Rz_Ry U_Lz_Ly U_Lz_Ry U_Rz_Ly U_Rz_Ry C_Lz C_Rz U_Lz U_Rz P).interface_f_x) = -((zInterfaceGradientExprs C_Lz_Lx
      C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz_Rx U_Lz_Lx U_Rz_Rx U_Rz_Lx C_Lz_Ly C_Lz_Ry C_Rz_Ly C_Rz_Ry U_Lz_Ry U_Lz_Ly U_Rz_Ry
      U_Rz_Ly C_Lz C_Rz U_Rz U_Lz P).interface_f_x) ∧
    ((zInterfaceGradientExprs C_Lz_Lx C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz_Lx U_Lz_Rx U_Rz_Lx U_Rz_Rx C_Lz_Ly C_Lz_Ry C_Rz_Ly
      C_Rz_Ry U_Lz_Ly U_Lz_Ry U_Rz_Ly U_Rz_Ry C_Lz C_Rz U_Lz U_Rz P).interface_f_y) = -((zInterfaceGradientExprs C_Lz_Lx
      C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz_Rx U_Lz_Lx U_Rz_Rx U_Rz_Lx C_Lz_Ly C_Lz_Ry C_Rz_Ly C_Rz_Ry U_Lz_Ry U_Lz_Ly U_Rz_Ry
      U_Rz_Ly C_Lz C_Rz U_Rz U_Lz P).interface_f_y) ∧
    ((zInterfaceGradientExprs C_Lz_Lx C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz_Lx U_Lz_Rx U_Rz_Lx U_Rz_Rx C_Lz_Ly C_Lz_Ry C_Rz_Ly
      C_Rz_Ry U_Lz_Ly U_Lz_Ry U_Rz_Ly U_Rz_Ry C_Lz C_Rz U_Lz U_Rz P).interface_f_z) = -((zInterfaceGradientExprs C_Lz_Lx
      C_Lz_Rx C_Rz_Lx C_Rz_Rx U_Lz_Rx U_Lz_Lx U_Rz_Rx U_Rz_Lx C_Lz_Ly C_Lz_Ry C_Rz_Ly C_Rz_Ry U_Lz_Ry U_Lz_Ly U_Rz_Ry
      U_Rz_Ly C_Lz C_Rz U_Rz U_Lz P).interface_f_z) := by
  constructor
  . simp [zInterfaceGradientExprs]
    ring_nf
  constructor
  . simp [zInterfaceGradientExprs]
    ring_nf
  . simp [zInterfaceGradientExprs]
    ring_nf

end advection_diffusion_iso_3d_centered_diffusion
                