import Mathlib

set_option linter.unusedVariables false
set_option linter.unusedTactic false
set_option maxHeartbeats 1000000

namespace linear_advection_2d

structure Coordinates where
  x : Real
  y : Real

structure State where
  f : Real

structure Parameters where
  ax : Real
  ay : Real

structure Flux where
  flux_f : Real

structure Fluxes2D where
  xFlux : Flux
  yFlux : Flux

noncomputable def xFluxExprs (C : Coordinates) (P : Parameters) (U : State) : Flux :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    flux_f := (ax * f)
  }

noncomputable def yFluxExprs (C : Coordinates) (P : Parameters) (U : State) : Flux :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    flux_f := (ay * f)
  }

noncomputable def fluxes (C : Coordinates) (P : Parameters) (U : State) : Fluxes2D :=
  {
    xFlux := xFluxExprs C P U
    yFlux := yFluxExprs C P U
  }

structure WaveSpeed where
  mu1 : Real

structure WaveSpeedEstimate2D where
  xWaveSpeed : WaveSpeed
  yWaveSpeed : WaveSpeed

noncomputable def xWaveSpeedExprs (C : Coordinates) (P : Parameters) (U : State) : WaveSpeed :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    mu1 := ax
  }

noncomputable def yWaveSpeedExprs (C : Coordinates) (P : Parameters) (U : State) : WaveSpeed :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    mu1 := ay
  }

noncomputable def wavespeeds (C : Coordinates) (P : Parameters) (U : State) : WaveSpeedEstimate2D :=
  {
    xWaveSpeed := xWaveSpeedExprs C P U
    yWaveSpeed := yWaveSpeedExprs C P U
  }

structure StateGradient where
  d_f : Real

structure FluxJacobian where
  flux_f : StateGradient

structure FluxJacobians2D where
  xFluxJacobian : FluxJacobian
  yFluxJacobian : FluxJacobian

noncomputable def xFluxJacobianExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobian :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
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
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    flux_f :=
      {
        d_f := ay
      }
  }

noncomputable def fluxJacobians (C : Coordinates) (P : Parameters) (U : State) : FluxJacobians2D :=
  {
    xFluxJacobian := xFluxJacobianExprs C P U
    yFluxJacobian := yFluxJacobianExprs C P U
  }

structure FluxJacobianEigen where
  lambda1 : Real

structure FluxJacobianEigendecomposition2D where
  xFluxJacobianEigen : FluxJacobianEigen
  yFluxJacobianEigen : FluxJacobianEigen

noncomputable def xFluxJacobianEigenExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigen :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    lambda1 := ax
  }

noncomputable def yFluxJacobianEigenExprs (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigen :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  {
    lambda1 := ay
  }

noncomputable def fluxJacobianEigenvalues (C : Coordinates) (P : Parameters) (U : State) : FluxJacobianEigendecomposition2D :=
  {
    xFluxJacobianEigen := xFluxJacobianEigenExprs C P U
    yFluxJacobianEigen := yFluxJacobianEigenExprs C P U
  }

theorem xHyperbolicity (C : Coordinates) (P : Parameters) (U : State) :
    (∃ r1 : Real, r1 = (xFluxJacobianEigenExprs C P U).lambda1) := by
  . refine ⟨(xFluxJacobianEigenExprs C P U).lambda1, rfl⟩

theorem yHyperbolicity (C : Coordinates) (P : Parameters) (U : State) :
    (∃ r1 : Real, r1 = (yFluxJacobianEigenExprs C P U).lambda1) := by
  . refine ⟨(yFluxJacobianEigenExprs C P U).lambda1, rfl⟩

theorem xWaveStability (C : Coordinates) (P : Parameters) (U : State) :
    ((abs (xWaveSpeedExprs C P U).mu1) ≥ (abs (xFluxJacobianEigenExprs C P U).lambda1)) := by
  . simp [xWaveSpeedExprs, xFluxJacobianEigenExprs]

theorem yWaveStability (C : Coordinates) (P : Parameters) (U : State) :
    ((abs (yWaveSpeedExprs C P U).mu1) ≥ (abs (yFluxJacobianEigenExprs C P U).lambda1)) := by
  . simp [yWaveSpeedExprs, yFluxJacobianEigenExprs]

structure SpatialGradient where
  f_x : Real
  f_y : Real

structure DiffusiveFlux where
  diffusive_flux_f : Real

structure DiffusiveFluxes2D where
  xDiffusiveFlux : DiffusiveFlux
  yDiffusiveFlux : DiffusiveFlux

noncomputable def xDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  {
    diffusive_flux_f := 0.0
  }

noncomputable def yDiffusiveFluxExprs (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFlux :=
  let x := C.x
  let y := C.y
  let ax := P.ax
  let ay := P.ay
  let f := U.f
  let f_x := DU.f_x
  let f_y := DU.f_y
  {
    diffusive_flux_f := 0.0
  }

noncomputable def diffusiveFluxes (C : Coordinates) (P : Parameters) (U : State) (DU : SpatialGradient) : DiffusiveFluxes2D :=
  {
    xDiffusiveFlux := xDiffusiveFluxExprs C P U DU
    yDiffusiveFlux := yDiffusiveFluxExprs C P U DU
  }

theorem xDiffusiveFluxConsistency (C : Coordinates) (U : State) (P : Parameters) :
    let DU : SpatialGradient :=
      {
        f_x := 0
        f_y := 0
      }
    ((xDiffusiveFluxExprs C P U DU).diffusive_flux_f) = 0 := by
  . simp only [xDiffusiveFluxExprs]
    norm_num

theorem yDiffusiveFluxConsistency (C : Coordinates) (U : State) (P : Parameters) :
    let DU : SpatialGradient :=
      {
        f_x := 0
        f_y := 0
      }
    ((yDiffusiveFluxExprs C P U DU).diffusive_flux_f) = 0 := by
  . simp only [yDiffusiveFluxExprs]
    norm_num

end linear_advection_2d

namespace linear_advection_2d_lax_friedrichs

open linear_advection_2d

structure Wave where
  wave_f : Real

structure WaveFamily where
  wave1 : Wave
  wave2 : Wave

structure Waves2D where
  xWaveFamily : WaveFamily
  yWaveFamily : WaveFamily

noncomputable def xWaveFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
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
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
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

noncomputable def waves (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Waves2D :=
  {
    xWaveFamily := xWaveFamilyExprs CL CR UL UR P
    yWaveFamily := yWaveFamilyExprs CL CR UL UR P
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

structure StateJump where
  jump_f : Real

structure StateJumps2D where
  xStateJump : StateJump
  yStateJump : StateJump

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

noncomputable def stateJumps (UL : State) (UR : State) : StateJumps2D :=
  {
    xStateJump := xStateJumpExprs UL UR
    yStateJump := yStateJumpExprs UL UR
  }

structure WaveSum where
  wave_sum_f : Real

structure WaveSums2D where
  xWaveSum : WaveSum
  yWaveSum : WaveSum

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

noncomputable def waveSums (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : WaveSums2D :=
  {
    xWaveSum := xWaveSumExprs CL CR UL UR P
    yWaveSum := yWaveSumExprs CL CR UL UR P
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

structure SpeedFamily where
  speed1 : Real
  speed2 : Real

structure Speeds2D where
  xSpeedFamily : SpeedFamily
  ySpeedFamily : SpeedFamily

noncomputable def xSpeedFamilyExprs (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : SpeedFamily :=
  let x_L := CL.x
  let x_R := CR.x
  let y_L := CL.y
  let y_R := CR.y
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
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
  let f_L := UL.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  {
    speed1 := -((abs ay))
    speed2 := (abs ay)
  }

noncomputable def speeds (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Speeds2D :=
  {
    xSpeedFamily := xSpeedFamilyExprs CL CR UL UR P
    ySpeedFamily := ySpeedFamilyExprs CL CR UL UR P
  }

structure Fluctuation where
  fluctuation_f : Real

structure Fluctuations2D where
  xFluctuation : Fluctuation
  yFluctuation : Fluctuation

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

noncomputable def leftFluctuations (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuations2D :=
  {
    xFluctuation := xLeftFluctuationExprs CL CR UL UR P
    yFluctuation := yLeftFluctuationExprs CL CR UL UR P
  }

noncomputable def rightFluctuations (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : Fluctuations2D :=
  {
    xFluctuation := xRightFluctuationExprs CL CR UL UR P
    yFluctuation := yRightFluctuationExprs CL CR UL UR P
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

structure FluxJump where
  flux_jump_f : Real

structure FluxJumps2D where
  xFluxJump : FluxJump
  yFluxJump : FluxJump

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

noncomputable def fluxJumps (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters): FluxJumps2D :=
  {
    xFluxJump := xFluxJumpExprs CL CR UL UR P
    yFluxJump := yFluxJumpExprs CL CR UL UR P
  }

structure FluctuationSum where
  fluctuation_sum_f : Real

structure FluctuationSums2D where
  xFluctuationSum : FluctuationSum
  yFluctuationSum : FluctuationSum

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

noncomputable def fluctuationSums (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters) : FluctuationSums2D :=
  {
    xFluctuationSum := xFluctuationSumExprs CL CR UL UR P
    yFluctuationSum := yFluctuationSumExprs CL CR UL UR P
  }

theorem xFluxConservative (CL : Coordinates) (CR : Coordinates)
    (UL : State) (UR : State) (P : Parameters)
    (Speed1Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      (-((abs ax)) ≠ 0)))
    (Speed2Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
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
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      (-((abs ay)) ≠ 0)))
    (Speed2Nondegenerate : (
      let x_L := CL.x
      let x_R := CR.x
      let y_L := CL.y
      let y_R := CR.y
      let f_L := UL.f
      let f_R := UR.f
      let ax := P.ax
      let ay := P.ay
      ((abs ay) ≠ 0))) : 
    ((yFluxJumpExprs CL CR UL UR P).flux_jump_f = (yFluctuationSumExprs CL CR UL UR P).fluctuation_sum_f) := by
  . simp [yFluxJumpExprs, yFluctuationSumExprs]
    simp [yFluxExprs, yLeftFluctuationExprs, yRightFluctuationExprs]
    ring_nf
    simp [yWaveFamilyExprs, ySpeedFamilyExprs]
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]
    ring_nf
    repeat field_simp [Speed1Nondegenerate, Speed2Nondegenerate]

end linear_advection_2d_lax_friedrichs

namespace linear_advection_2d_minmod

open linear_advection_2d

structure Reconstruction where
  reconstruction_f : Real

structure Reconstruction2D where
  xReconstruction : Reconstruction
  yReconstruction : Reconstruction

noncomputable def xLeftReconstructionExprs (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction :=
  let x_L := CL.x
  let x := C.x
  let x_R := CR.x
  let y_L := CL.y
  let y := C.y
  let y_R := CR.y
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
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
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
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
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
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
  let f_L := UL.f
  let f := U.f
  let f_R := UR.f
  let ax := P.ax
  let ay := P.ay
  {
    reconstruction_f := (f + ((1 / 2 : Real) * ((max 0.0 (min (f - f_L) (f_R - f))) + (min 0.0 (max (f - f_L) (f_R -
      f))))))
  }

noncomputable def leftReconstruction (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction2D :=
  {
    xReconstruction := xLeftReconstructionExprs CL C CR UL U UR P
    yReconstruction := yLeftReconstructionExprs CL C CR UL U UR P
  }

noncomputable def rightReconstruction (CL : Coordinates) (C : Coordinates) (CR : Coordinates)
    (UL : State) (U : State) (UR : State) (P : Parameters) : Reconstruction2D :=
  {
    xReconstruction := xRightReconstructionExprs CL C CR UL U UR P
    yReconstruction := yRightReconstructionExprs CL C CR UL U UR P
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

end linear_advection_2d_minmod
                