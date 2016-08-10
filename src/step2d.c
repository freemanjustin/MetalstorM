#include "MetalstorM.h"

void step2d(e *E){

  int i,j;

  ptsk=3-kstp
  CORRECTOR_2D_STEP=.not.PREDICTOR_2D_STEP(ng)
  //
  //-----------------------------------------------------------------------
  //  Compute total depth (m) and vertically integrated mass fluxes.
  //-----------------------------------------------------------------------
  //
  for(j=JstrVm2-1; j<Jendp2; j++){
    for(i=IstrUm2-1; i<Iendp2; i++){
      Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
    }
  }

  for(j=JstrVm2-1; j<Jendp2; j++){
    for(i=IstrUm2; i<Iendp2; i++){
      cff=0.5_r8*on_u(i,j)
      cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
      DUon(i,j)=ubar(i,j,krhs)*cff1
    }
  }
  for(j=JstrVm2; j<Jendp2; j++){
    for(i=IstrUm2-1; i<Iendp2; i++){
      cff=0.5_r8*om_v(i,j)
      cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
      DVom(i,j)=vbar(i,j,krhs)*cff1
    }
  }
  //
  //  Set vertically integrated mass fluxes DUon and DVom along the open
  //  boundaries in such a way that the integral volume is conserved.
  //
        if (ANY(VolCons(:,ng))) {
          CALL set_DUV_bc_tile (ng, tile,                                 &
       &                        LBi, UBi, LBj, UBj,                       &
       &                        IminS, ImaxS, JminS, JmaxS,               &
       &                        krhs,                                     &
       &                        umask, vmask,                             &
       &                        om_v, on_u,                               &
       &                        ubar, vbar,                               &
       &                        Drhs, DUon, DVom)
        }
  //
  //-----------------------------------------------------------------------
  //  Compute new wet/dry masks.
  //-----------------------------------------------------------------------
  //
        CALL wetdry_tile (ng, tile,                                       &
       &                  LBi, UBi, LBj, UBj,                             &
       &                  IminS, ImaxS, JminS, JmaxS,                     &
       &                  pmask, rmask, umask, vmask,                     &
       &                  h, zeta(:,:,kstp),                              &
       &                  pmask_wet, pmask_full,                          &
       &                  rmask_wet, rmask_full,                          &
       &                  umask_wet, umask_full,                          &
       &                  vmask_wet, vmask_full)
  //
  //  for(not perform the actual time stepping during the auxiliary
  //  (nfast(ng)+1) time step.
  //
        if (iif(ng).gt.nfast(ng)) RETURN
  //
  //=======================================================================
  //  Time step free-surface equation.
  //=======================================================================
  //
  //  During the first time-step, the predictor step is Forward-Euler
  //  and the corrector step is Backward-Euler. Otherwise, the predictor
  //  step is Leap-frog and the corrector step is Adams-Moulton.
  //
        if (iic(ng).eq.ntfirst(ng)){
          cff1=dtfast(ng)
          for(j=JstrV-1; j<Jend; j++){
            for(i=IstrU-1; i<Iend; i++){
              rhs_zeta(i,j)=(DUon(i,j)-DUon(i+1,j))+(DVom(i,j)-DVom(i,j+1))
              zeta_new(i,j)=zeta(i,j,kstp)+pm(i,j)*pn(i,j)*cff1*rhs_zeta(i,j)
              zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
              Dnew(i,j)=zeta_new(i,j)+h(i,j)

              zwrk(i,j)=0.5_r8*(zeta(i,j,kstp)+zeta_new(i,j))
              gzeta(i,j)=zwrk(i,j)
              gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
            }
          }
        }
        else if (PREDICTOR_2D_STEP(ng)){
          cff1=2.0_r8*dtfast(ng)
          cff4=4.0_r8/25.0_r8
          cff5=1.0_r8-2.0_r8*cff4
          for(j=JstrV-1; j<Jend; j++){
            for(i=IstrU-1; i<Iend; i++){
              rhs_zeta(i,j)=(DUon(i,j)-DUon(i+1,j))+(DVom(i,j)-DVom(i,j+1))
              zeta_new(i,j)=zeta(i,j,kstp)+pm(i,j)*pn(i,j)*cff1*rhs_zeta(i,j)
              zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
              Dnew(i,j)=zeta_new(i,j)+h(i,j)
  //
              zwrk(i,j)=cff5*zeta(i,j,krhs)+cff4*(zeta(i,j,kstp)+zeta_new(i,j))
              gzeta(i,j)=zwrk(i,j)
              gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
            }
          }
        }
        else if(CORRECTOR_2D_STEP){
          cff1=dtfast(ng)*5.0_r8/12.0_r8
          cff2=dtfast(ng)*8.0_r8/12.0_r8
          cff3=dtfast(ng)*1.0_r8/12.0_r8
          cff4=2.0_r8/5.0_r8
          cff5=1.0_r8-cff4
          for(j=JstrV-1; j<Jend; j++){
            for(i=IstrU-1; i<Iend; i++){
              cff=cff1*((DUon(i,j)-DUon(i+1,j))+ (DVom(i,j)-DVom(i,j+1)))
              zeta_new(i,j)=zeta(i,j,kstp)+pm(i,j)*pn(i,j)*
                            (cff+cff2*rzeta(i,j,kstp)-cff3*rzeta(i,j,ptsk))
              zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
              Dnew(i,j)=zeta_new(i,j)+h(i,j)
  //
              zwrk(i,j)=cff5*zeta_new(i,j)+cff4*zeta(i,j,krhs)
              gzeta(i,j)=zwrk(i,j)
              gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
            }
          }
        }
  //
  //  Load new free-surface values into shared array at both predictor
  //  and corrector steps.
  //  Modify new free-surface to Ensure that depth is > Dcrit for masked
  //  cells.
  //
        for(j=Jstr; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            zeta(i,j,knew)=zeta_new(i,j)
            zeta(i,j,knew)=zeta(i,j,knew)+(Dcrit(ng)-h(i,j))*(1.0_r8-rmask(i,j))
          }
        }
  //
  //  if predictor step, load right-side-term into shared array.
  //
        if (PREDICTOR_2D_STEP(ng)){
          for(j=Jstr; j<Jend; j++){
            for(i=Istr; i<Iend; i++){
              rzeta(i,j,krhs)=rhs_zeta(i,j)
            }
          }
          if (EWperiodic(ng).or.NSperiodic(ng)){
            CALL exchange_r2d_tile (ng, tile,                             &
       &                            LBi, UBi, LBj, UBj,                   &
       &                            rzeta(:,:,krhs))
          }
      }
  //
  //  Apply mass point sources (volume vertical influx), if any.
  //
        if (LwSrc(ng)){
          for(is=1; is<Nsrc(ng); is++){
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            if (((IstrR.le.i).and.(i.le.IendR)) && ((JstrR.le.j).and.(j.le.JendR))){
              zeta(i,j,knew)=zeta(i,j,knew)+                              &
       &                     SOURCES(ng)%Qbar(is)*                        &
       &                     pm(i,j)*pn(i,j)*dtfast(ng)
            }
          }
       }
  //
  //  Set free-surface lateral boundary conditions.
  //
        CALL zetabc_tile (ng, tile,                                       &
       &                  LBi, UBi, LBj, UBj,                             &
       &                  IminS, ImaxS, JminS, JmaxS,                     &
       &                  krhs, kstp, knew,                               &
       &                  zeta)
        if (EWperiodic(ng).or.NSperiodic(ng)){
          CALL exchange_r2d_tile (ng, tile,                               &
       &                          LBi, UBi, LBj, UBj,                     &
       &                          zeta(:,:,knew))
        }
  //
  //=======================================================================
  //  Compute right-hand-side for the 2D momentum equations.
  //=======================================================================
  //
  //-----------------------------------------------------------------------
  //  Compute pressure gradient terms.
  //-----------------------------------------------------------------------
  //
        cff1=0.5_r8*g
        cff2=1.0_r8/3.0_r8
        fac3=0.5_r8*100.0_r8/rho0
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            rhs_ubar(i,j)=cff1*on_u(i,j)* ((h(i-1,j)+ h(i ,j))*
                          (gzeta(i-1,j)-gzeta(i  ,j))+(gzeta2(i-1,j)-gzeta2(i,j)))
            rhs_ubar(i,j)=rhs_ubar(i,j)+fac3*on_u(i,j)*(h(i-1,j)+h(i,j)
                          +gzeta(i-1,j)+gzeta(i,j))
                          *(Pair(i-1,j)-Pair(i,j))
          }
          if (j>=JstrV){
            for(i=Istr; i<Iend; i++){
              rhs_vbar(i,j)=cff1*om_v(i,j)*((h(i,j-1)+h(i,j  ))
                            *(gzeta(i,j-1)-gzeta(i,j))+ (gzeta2(i,j-1)-gzeta2(i,j  )))
              rhs_vbar(i,j)=rhs_vbar(i,j)+fac3*om_v(i,j)*(h(i,j-1)+h(i,j)
                            +gzeta(i,j-1)+gzeta(i,j))
                            *(Pair(i,j-1)-Pair(i,j))
            }
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Add in horizontal advection of momentum.
  //-----------------------------------------------------------------------
  //
  //  Fourth-order, centered differences advection.
  //
        for(j=Jstr; j<Jend; j++){
          for(i=IstrUm1; i<Iendp1; i++){
            grad (i,j)=ubar(i-1,j,krhs)-2.0_r8*ubar(i,j,krhs)+ubar(i+1,j,krhs)
            Dgrad(i,j)=DUon(i-1,j)-2.0_r8*DUon(i,j)+DUon(i+1,j)
          }
        }
        if (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))){
          if (DOMAIN(ng)%Western_Edge(tile)){
            for(j=Jstr; j<Jend; j++){
              grad (Istr,j)=grad (Istr+1,j)
              Dgrad(Istr,j)=Dgrad(Istr+1,j)
            }
          }
        }
        if (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))){
          if (DOMAIN(ng)%Eastern_Edge(tile)){
            for(j=Jstr; j<Jend; j++){
              grad (Iend+1,j)=grad (Iend,j)
              Dgrad(Iend+1,j)=Dgrad(Iend,j)
            }
          }
        }
        cff=1.0_r8/6.0_r8
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU-1; i<Iend; i++){
            UFx(i,j)=0.25_r8*(ubar(i  ,j,krhs)+ubar(i+1,j,krhs)
                      -cff*(grad (i,j)+grad (i+1,j)))*(DUon(i,j)+DUon(i+1,j)
                      -cff*(Dgrad(i,j)+Dgrad(i+1,j)))
          }
        }
  //
        for(j=Jstrm1; j<Jendp1; j++){
          for(i=IstrU; i<Iend; i++){
            grad(i,j)=ubar(i,j-1,krhs)-2.0_r8*ubar(i,j,krhs)+ubar(i,j+1,krhs)
          }
        }
        if (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))){
          if (DOMAIN(ng)%Southern_Edge(tile)){
            for(i=IstrU; i<Iend; i++){
              grad(i,Jstr-1)=grad(i,Jstr)
            }
          }
        }
        if (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) {
          if (DOMAIN(ng)%Northern_Edge(tile)) {
            for(i=IstrU; i<Iend; i++){
              grad(i,Jend+1)=grad(i,Jend)
            }
          }
        }
        for(j=Jstr; j<Jend+1; j++){
          for(i=IstrU-1; i<Iend; i++){
            Dgrad(i,j)=DVom(i-1,j)-2.0_r8*DVom(i,j)+DVom(i+1,j)
          }
        }
        cff=1.0_r8/6.0_r8
        for(j=Jstr; j<Jend+1; j++){
          for(i=IstrU; i<Iend; i++){
            UFe(i,j)=0.25_r8*(ubar(i,j  ,krhs)+
                      ubar(i,j-1,krhs)-
                      cff*(grad (i,j)+grad (i,j-1)))*
                      (DVom(i,j)+DVom(i-1,j)-
                      cff*(Dgrad(i,j)+Dgrad(i-1,j)))
          }
        }
  //
        for(j=JstrV; j<Jend; j++){
          for(i=Istrm1; i<Iendp1; i++){
            grad(i,j)=vbar(i-1,j,krhs)-2.0_r8*vbar(i,j,krhs)+
                      vbar(i+1,j,krhs)
          }
        }
        if (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) {
          if (DOMAIN(ng)%Western_Edge(tile)) {
            for(j=JstrV; j<Jend; j++){
              grad(Istr-1,j)=grad(Istr,j)
            }
          }
        }
        if (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) {
          if (DOMAIN(ng)%Eastern_Edge(tile)) {
            for(j=JstrV; j<Jend; j++){
              grad(Iend+1,j)=grad(Iend,j)
            }
          }
        }
        for(j=JstrV-1; j<Jend; j++){
          for(i=Istr; i<Iend+1; i++){
            Dgrad(i,j)=DUon(i,j-1)-2.0_r8*DUon(i,j)+DUon(i,j+1)
          }
        }
        cff=1.0_r8/6.0_r8
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend+1; i++){
            VFx(i,j)=0.25_r8*(vbar(i  ,j,krhs)+
                      vbar(i-1,j,krhs)-
                      cff*(grad (i,j)+grad (i-1,j)))*
                      (DUon(i,j)+DUon(i,j-1)-
                      cff*(Dgrad(i,j)+Dgrad(i,j-1)))
          }
        }
  //
        for(j=JstrVm1; j<Jendp1; j++){
          for(i=Istr; i<Iend; i++){
            grad(i,j)=vbar(i,j-1,krhs)-2.0_r8*vbar(i,j,krhs)+
                      vbar(i,j+1,krhs)
            Dgrad(i,j)=DVom(i,j-1)-2.0_r8*DVom(i,j)+DVom(i,j+1)
          }
        }
        if (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) {
          if (DOMAIN(ng)%Southern_Edge(tile)) {
            for(i=Istr; i<Iend; i++){
              grad (i,Jstr)=grad (i,Jstr+1)
              Dgrad(i,Jstr)=Dgrad(i,Jstr+1)
            }
          }
        }
        if (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) {
          if (DOMAIN(ng)%Northern_Edge(tile)) {
            for(i=Istr; i<Iend; i++){
              grad (i,Jend+1)=grad (i,Jend)
              Dgrad(i,Jend+1)=Dgrad(i,Jend)
            }
          }
        }
        cff=1.0_r8/6.0_r8
        for(j=JstrV-1; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            VFe(i,j)=0.25_r8*(vbar(i,j  ,krhs)+
                          vbar(i,j+1,krhs)-
                          cff*(grad (i,j)+grad (i,j+1)))*
                          (DVom(i,j)+DVom(i,j+1)-
                          cff*(Dgrad(i,j)+Dgrad(i,j+1)))
          }
        }
  //
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=UFe(i,j+1)-UFe(i,j)
            fac=cff1+cff2
            rhs_ubar(i,j)=rhs_ubar(i,j)-fac
          }
        }
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            cff1=VFx(i+1,j)-VFx(i,j)
            cff2=VFe(i,j)-VFe(i,j-1)
            fac=cff1+cff2
            rhs_vbar(i,j)=rhs_vbar(i,j)-fac
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Add in Coriolis term.
  //-----------------------------------------------------------------------
  //
        for(j=JstrV-1; j<Jend; j++){
          for(i=IstrU-1; i<Iend; i++){
            cff=0.5_r8*Drhs(i,j)*fomn(i,j)
            UFx(i,j)=cff*(vbar(i,j  ,krhs)+vbar(i,j+1,krhs))
            VFe(i,j)=cff*(ubar(i  ,j,krhs)+ubar(i+1,j,krhs))
          }
        }
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            rhs_ubar(i,j)=rhs_ubar(i,j)+fac1
          }
        }
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            fac1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rhs_vbar(i,j)=rhs_vbar(i,j)-fac1
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Add in curvilinear transformation terms.
  //-----------------------------------------------------------------------
  //
        for(j=JstrV-1; j<Jend; j++){
          for(i=IstrU-1; i<Iend; i++){
            cff1=0.5_r8*(vbar(i,j  ,krhs)+vbar(i,j+1,krhs))
            cff2=0.5_r8*(ubar(i  ,j,krhs)+ubar(i+1,j,krhs))
            cff3=cff1*dndx(i,j)
            cff4=cff2*dmde(i,j)
            cff=Drhs(i,j)*(cff3-cff4)
            UFx(i,j)=cff*cff1
            VFe(i,j)=cff*cff2
          }
        }
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            rhs_ubar(i,j)=rhs_ubar(i,j)+fac1
          }
        }
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            fac1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rhs_vbar(i,j)=rhs_vbar(i,j)-fac1
          }
        }
  //
  //-----------------------------------------------------------------------
  //  if horizontal mixing, compute total depth at PSI-points.
  //-----------------------------------------------------------------------
  //
        for(j=Jstr; j<Jend+1; j++){
          for(i=Istr; i<Iend+1; i++){
            Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+Drhs(i,j-1)+Drhs(i-1,j-1))
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Add in horizontal harmonic viscosity.
  //-----------------------------------------------------------------------
  //
  //  Compute flux-components of the horizontal divergence of the stress
  //  tensor (m5/s2) in XI- and ETA-directions.
  //
        for(j=JstrV-1; j<Jend; j++){
          for(i=IstrU-1; i<Iend; i++){
            cff=visc2_r(i,j)*Drhs(i,j)*0.5_r8*
                  (pmon_r(i,j)*
                  ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-
                  (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-
                  pnom_r(i,j)*
                  ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-
                  (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))
            UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
            VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
          }
        }
        for(j=Jstr; j<Jend+1; j++){
          for(i=Istr; i<Iend+1; i++){
            cff=visc2_p(i,j)*Drhs_p(i,j)*0.5_r8*
                  (pmon_p(i,j)*
                  ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-
                  (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+
                  pnom_p(i,j)*
                  ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-
                  (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))
            cff=cff*pmask(i,j)
            cff=cff*pmask_wet(i,j)
            UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
            VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
          }
        }
  //
  //  Add in harmonic viscosity.
  //
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))
            cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j))
            fac=cff1+cff2
            rhs_ubar(i,j)=rhs_ubar(i,j)+fac
          }
        }
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
            cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
            fac=cff1-cff2
            rhs_vbar(i,j)=rhs_vbar(i,j)+fac
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Add in bottom stress.
  //-----------------------------------------------------------------------
  //
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            fac=bustr(i,j)*om_u(i,j)*on_u(i,j)
            rhs_ubar(i,j)=rhs_ubar(i,j)-fac
          }
        }
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            fac=bvstr(i,j)*om_v(i,j)*on_v(i,j)
            rhs_vbar(i,j)=rhs_vbar(i,j)-fac
          }
        }

  //
  //-----------------------------------------------------------------------
  //  Add in surface momentum stress.
  //-----------------------------------------------------------------------
  //
        for(j=Jstr; j<Jend; j++){
          for(i=IstrU; i<Iend; i++){
            fac=sustr(i,j)*om_u(i,j)*on_u(i,j)
            rhs_ubar(i,j)=rhs_ubar(i,j)+fac
          }
        }
        for(j=JstrV; j<Jend; j++){
          for(i=Istr; i<Iend; i++){
            fac=svstr(i,j)*om_v(i,j)*on_v(i,j)
            rhs_vbar(i,j)=rhs_vbar(i,j)+fac
          }
        }
  //
  //=======================================================================
  //  Time step 2D momentum equations.
  //=======================================================================
  //
  //  Compute total water column depth.
  //
        for(j=JstrV-1; j<Jend; j++){
          for(i=IstrU-1; i<Iend; i++){
            Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
          }
        }
  //
  //  During the first time-step, the predictor step is Forward-Euler
  //  and the corrector step is Backward-Euler. Otherwise, the predictor
  //  step is Leap-frog and the corrector step is Adams-Moulton.
  //
        if (iic(ng).eq.ntfirst(ng)) {
          cff1=0.5_r8*dtfast(ng)
          cff2=1.0_r8/cff1
          for(j=Jstr; j<Jend; j++){
            for(i=IstrU; i<Iend; i++){
              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
              ubar(i,j,knew)=(ubar(i,j,kstp)*
                            (Dstp(i,j)+Dstp(i-1,j))+
                            cff*cff1*rhs_ubar(i,j))*fac
              ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
              cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
              cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
              cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
              ubar(i,j,knew)=ubar(i,j,knew)*cff7
              rhs_ubar(i,j)=rhs_ubar(i,j)*cff7
            }
          }
          for(j=JstrV; j<Jend; j++){
            for(i=Istr; i<Iend; i++){
              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
              vbar(i,j,knew)=(vbar(i,j,kstp)*
                              (Dstp(i,j)+Dstp(i,j-1))+
                              cff*cff1*rhs_vbar(i,j))*fac
              vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
              cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
              cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
              cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
              vbar(i,j,knew)=vbar(i,j,knew)*cff7
              rhs_vbar(i,j)=rhs_vbar(i,j)*cff7
            }
          }
        }
        else if (PREDICTOR_2D_STEP(ng)) {
          cff1=dtfast(ng)
          cff2=1.0_r8/cff1
          for(j=Jstr; j<Jend; j++){
            for(i=IstrU; i<Iend; i++){
              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
              ubar(i,j,knew)=(ubar(i,j,kstp)*
                              (Dstp(i,j)+Dstp(i-1,j))+
                              cff*cff1*rhs_ubar(i,j))*fac
              ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
              cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
              cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
              cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
              ubar(i,j,knew)=ubar(i,j,knew)*cff7
              rhs_ubar(i,j)=rhs_ubar(i,j)*cff7
            }
          }
          for(j=JstrV; j<Jend; j++){
            for(i=Istr; i<Iend; i++){
              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
              vbar(i,j,knew)=(vbar(i,j,kstp)*
                              (Dstp(i,j)+Dstp(i,j-1))+
                              cff*cff1*rhs_vbar(i,j))*fac
              vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
              cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
              cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
              cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
              vbar(i,j,knew)=vbar(i,j,knew)*cff7
              rhs_vbar(i,j)=rhs_vbar(i,j)*cff7
            }
          }
        }
        else if (CORRECTOR_2D_STEP) {
          cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
          cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
          cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
          cff4=1.0_r8/cff1
          for(j=Jstr; j<Jend; j++){
            for(i=IstrU; i<Iend; i++){
              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
              ubar(i,j,knew)=(ubar(i,j,kstp)*
                              (Dstp(i,j)+Dstp(i-1,j))+
                              cff*(cff1*rhs_ubar(i,j)+
                              cff2*rubar(i,j,kstp)-
                              cff3*rubar(i,j,ptsk)))*fac
              ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
              cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
              cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
              cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
              ubar(i,j,knew)=ubar(i,j,knew)*cff7
              rhs_ubar(i,j)=rhs_ubar(i,j)*cff7
            }
          }
          for(j=JstrV; j<Jend; j++){
            for(i=Istr; i<Iend; i++){
              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
              vbar(i,j,knew)=(vbar(i,j,kstp)*
                              (Dstp(i,j)+Dstp(i,j-1))+
                              cff*(cff1*rhs_vbar(i,j)+
                              cff2*rvbar(i,j,kstp)-
                              cff3*rvbar(i,j,ptsk)))*fac
              vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
              cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
              cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
              cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
              vbar(i,j,knew)=vbar(i,j,knew)*cff7
              rhs_vbar(i,j)=rhs_vbar(i,j)*cff7
            }
          }
        }
  //
  //  if predictor step, load right-side-term into shared arrays for
  //  future use during the subsequent corrector step.
  //
        if (PREDICTOR_2D_STEP(ng)) {
          for(j=Jstr; j<Jend; j++){
            for(i=IstrU; i<Iend; i++){
              rubar(i,j,krhs)=rhs_ubar(i,j)
            }
          }
          for(j=JstrV; j<Jend; j++){
            for(i=Istr; i<Iend; i++){
              rvbar(i,j,krhs)=rhs_vbar(i,j)
            }
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Apply lateral boundary conditions.
  //-----------------------------------------------------------------------
  //
        CALL u2dbc_tile (ng, tile,                                        &
       &                 LBi, UBi, LBj, UBj,                              &
       &                 IminS, ImaxS, JminS, JmaxS,                      &
       &                 krhs, kstp, knew,                                &
       &                 ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                        &
       &                 LBi, UBi, LBj, UBj,                              &
       &                 IminS, ImaxS, JminS, JmaxS,                      &
       &                 krhs, kstp, knew,                                &
       &                 ubar, vbar, zeta)
  //
  //  Compute integral mass flux across open boundaries and adjust
  //  for volume conservation.
  //
        if (ANY(VolCons(:,ng))) {
          CALL obc_flux_tile (ng, tile,                                   &
       &                      LBi, UBi, LBj, UBj,                         &
       &                      IminS, ImaxS, JminS, JmaxS,                 &
       &                      knew,                                       &
       &                      umask, vmask,                               &
       &                      h, om_v, on_u,                              &
       &                      ubar, vbar, zeta)
        }
  //
  //-----------------------------------------------------------------------
  //  Apply momentum transport point sources (like river runoff), if any.
  //-----------------------------------------------------------------------
  //
        if (LuvSrc(ng)) {
          for(is=1; is<Nsrc(ng); is++){
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            if (((IstrR.le.i).and.(i.le.IendR))&&((JstrR.le.j).and.(j.le.JendR))) {
              if (INT(SOURCES(ng)%Dsrc(is))==0) {
                cff=1.0_r8/(on_u(i,j)*
                    0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+
                    zeta(i  ,j,knew)+h(i  ,j)))
                ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
              }
              else{
                cff=1.0_r8/(om_v(i,j)*
                    0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+
                    zeta(i,j  ,knew)+h(i,j  )))
                vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
              }
            }
          }
        }
  //
  //-----------------------------------------------------------------------
  //  Exchange boundary information.
  //-----------------------------------------------------------------------
  //
        if (EWperiodic(ng) || NSperiodic(ng)) {
          CALL exchange_u2d_tile (ng, tile,                               &
       &                          LBi, UBi, LBj, UBj,                     &
       &                          ubar(:,:,knew))
          CALL exchange_v2d_tile (ng, tile,                               &
       &                          LBi, UBi, LBj, UBj,                     &
       &                          vbar(:,:,knew))
        }

}
