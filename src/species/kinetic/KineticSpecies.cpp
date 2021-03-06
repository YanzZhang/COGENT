#include "FourthOrderUtil.H"
#include "KineticSpecies.H"
#include "MomentOp.H"
#include "Misc.H"
#include "Directions.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluxSurface.H"
#include "FourthOrderUtil.H"
#include "MagGeom.H"
#include "MagBlockCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

KineticSpecies::KineticSpecies(
         const string&                   a_name,
         const Real                      a_mass,
         const Real                      a_charge,
         const RefCountedPtr<PhaseGeom>& a_geometry
   )
  : m_geometry( a_geometry ),
    m_name( a_name ),
    m_mass( a_mass ),
    m_charge( a_charge ),
    m_moment_op( MomentOp::instance() )
{
}


KineticSpecies::KineticSpecies( const KineticSpecies& a_foo )
  : m_geometry( a_foo.m_geometry ),
    m_name( a_foo.m_name ),
    m_mass( a_foo.m_mass ),
    m_charge( a_foo.m_charge ),
    m_moment_op( MomentOp::instance() )
{
   m_dist_func.define( a_foo.m_dist_func );
}

void KineticSpecies::numberDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   m_moment_op.compute( a_rho, *this, DensityKernel() );
}


void KineticSpecies::massDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   m_moment_op.compute( a_rho, *this, MassDensityKernel() );
}


void KineticSpecies::chargeDensity( CFG::LevelData<CFG::FArrayBox>& a_rho ) const
{
   m_moment_op.compute( a_rho, *this, ChargeDensityKernel() );
}


void KineticSpecies::momentumDensity( CFG::LevelData<CFG::FArrayBox>& a_momentum,
				      const LevelData<FluxBox>& a_field,
				      const double larmor  ) const
{
   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();

   CFG::LevelData<CFG::FArrayBox> guidingCenterMom;
   guidingCenterMom.define(a_momentum); 
   m_moment_op.compute( guidingCenterMom, *this, MomentumDensityKernel(a_field) );
  

   CFG::LevelData<CFG::FArrayBox> magnetization(a_momentum.disjointBoxLayout(), 1, a_momentum.ghostVect());
   m_moment_op.compute( magnetization, *this, MagnetizationKernel() );

   CFG::LevelData<CFG::FArrayBox> magnetization_grown(a_momentum.disjointBoxLayout(), 1, a_momentum.ghostVect()+2*CFG::IntVect::Unit);
   
   //Extrapolating magnetization moment into the two layers of ghost cells
   for (CFG::DataIterator dit(magnetization.dataIterator()); dit.ok(); ++dit) {
     int block_number( coords.whichBlock( mag_geom.gridsFull()[dit] ) );
     const CFG::MagBlockCoordSys& block_coord_sys = (const CFG::MagBlockCoordSys&)(*(coords.getCoordSys(block_number)));
     const CFG::ProblemDomain& block_domain = block_coord_sys.domain();
     magnetization_grown[dit].copy(magnetization[dit]);
     fourthOrderCellExtrapAtDomainBdry(magnetization_grown[dit], block_domain, mag_geom.gridsFull()[dit]);
   }
   mag_geom.fillInternalGhosts(magnetization_grown);
   
   //Compute radial component of the magnetization gradient (part of curlM calculation)
   CFG::LevelData<CFG::FArrayBox> gradM_mapped(a_momentum.disjointBoxLayout(), CFG_DIM,
                                               a_momentum.ghostVect()+CFG::IntVect::Unit);
  
   mag_geom.computeMappedPoloidalGradientWithGhosts(magnetization_grown,gradM_mapped,2);
   
   CFG::LevelData<CFG::FArrayBox> gradM_phys;
   gradM_phys.define(gradM_mapped);
   mag_geom.unmapPoloidalGradient(gradM_mapped, gradM_phys);
   
   //Get magnetic geometry data
   const CFG::LevelData<CFG::FArrayBox>& bunit = mag_geom.getCCBFieldDir();
   const CFG::LevelData<CFG::FArrayBox>& curlb = mag_geom.getCCCurlBFieldDir();

   double fac = larmor / 2.0 / m_charge;

   CFG::LevelData<CFG::FArrayBox> Mcurlb;
   Mcurlb.define(a_momentum);
   CFG::LevelData<CFG::FArrayBox> bxgradM;
   bxgradM.define(a_momentum);

   for (CFG::DataIterator dit(a_momentum.dataIterator()); dit.ok(); ++dit) {
      
     //Compute M*curl(b) term
     Mcurlb[dit].copy(magnetization[dit],0,0,1);
     Mcurlb[dit].copy(magnetization[dit],0,1,1);
     Mcurlb[dit].mult(curlb[dit],0,0,1);
     Mcurlb[dit].mult(curlb[dit],2,1,1);
     Mcurlb[dit].mult(fac);

     //Compute b x gradM term
     bxgradM[dit].copy(gradM_phys[dit],1,0,1);
     bxgradM[dit].mult(bunit[dit],1,0,1);
     bxgradM[dit].mult(-1.0,0,1);

     bxgradM[dit].copy(gradM_phys[dit],0,1,1);
     bxgradM[dit].mult(bunit[dit],1,1,1);

     bxgradM[dit].mult(fac);
 
     //Add the poloidal component of the guiding center velocity and the curl of magnetization
     a_momentum[dit].copy(guidingCenterMom[dit],0,0,CFG_DIM);
     a_momentum[dit].plus(Mcurlb[dit],0,0,CFG_DIM);
     a_momentum[dit].plus(bxgradM[dit],0,0,CFG_DIM);
   }



}

void KineticSpecies::ParallelMomentum( CFG::LevelData<CFG::FArrayBox>& a_ParallelMom ) const
{
   m_moment_op.compute( a_ParallelMom, *this, ParallelMomKernel() );
}

void KineticSpecies::PoloidalMomentum( CFG::LevelData<CFG::FArrayBox>& a_PoloidalMom,
                                       const LevelData<FluxBox>& field,
                                       const double larmor  ) const 
{


   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();

   CFG::LevelData<CFG::FArrayBox> guidingCenterPoloidalMom;
   guidingCenterPoloidalMom.define(a_PoloidalMom); 
   m_moment_op.compute( guidingCenterPoloidalMom, *this, GuidingCenterPoloidalMomKernel(field) );
   
   CFG::LevelData<CFG::FArrayBox> magnetization;
   magnetization.define(a_PoloidalMom);
   m_moment_op.compute( magnetization, *this, MagnetizationKernel() );

   CFG::LevelData<CFG::FArrayBox> magnetization_grown(a_PoloidalMom.disjointBoxLayout(),1,
                                                      a_PoloidalMom.ghostVect()+2*CFG::IntVect::Unit);
   
   //Extrapolating magnetization moment into the two layers of ghost cells
   for (CFG::DataIterator dit(magnetization.dataIterator()); dit.ok(); ++dit) {
      int block_number( coords.whichBlock( mag_geom.gridsFull()[dit] ) );
      const CFG::MagBlockCoordSys& block_coord_sys = (const CFG::MagBlockCoordSys&)(*(coords.getCoordSys(block_number)));
      const CFG::ProblemDomain& block_domain = block_coord_sys.domain();
      magnetization_grown[dit].copy(magnetization[dit]);
      fourthOrderCellExtrapAtDomainBdry(magnetization_grown[dit], block_domain, mag_geom.gridsFull()[dit]);
   }
   mag_geom.fillInternalGhosts(magnetization_grown);
   
   //Compute radial component of the magnetization gradient (part of curlM calculation)
   CFG::LevelData<CFG::FArrayBox> gradM_mapped(a_PoloidalMom.disjointBoxLayout(), 2,
                                               a_PoloidalMom.ghostVect()+CFG::IntVect::Unit);
  
   mag_geom.computeMappedPoloidalGradientWithGhosts(magnetization_grown,gradM_mapped,2);
   
   CFG::LevelData<CFG::FArrayBox> gradM_phys;
   gradM_phys.define(gradM_mapped);
   mag_geom.unmapPoloidalGradient(gradM_mapped, gradM_phys);
   
   CFG::LevelData<CFG::FArrayBox> gradM_r(gradM_phys.disjointBoxLayout(),1, gradM_phys.ghostVect());
   mag_geom.computeRadialProjection(gradM_r, gradM_phys);
   
   //Get magnetic geometry data
   const CFG::LevelData<CFG::FArrayBox>& bunit = mag_geom.getCCBFieldDir();
   const CFG::LevelData<CFG::FArrayBox>& curlb = mag_geom.getCCCurlBFieldDir();
   CFG::LevelData<CFG::FArrayBox> curlb_pol(curlb.disjointBoxLayout(),1,curlb.ghostVect());
   mag_geom.computePoloidalProjection(curlb_pol, curlb);

   
   double fac = larmor / 2.0 / m_charge;

   for (CFG::DataIterator dit(a_PoloidalMom.dataIterator()); dit.ok(); ++dit) {
      
      //Compute curl of Magnetization vector
      curlb_pol[dit].mult(magnetization[dit]);
      curlb_pol[dit].mult(fac);
      gradM_r[dit].mult(bunit[dit],1,0,1);
      gradM_r[dit].mult(fac);
      
      //Add the poloidal component of the guiding center velocity and the curl of magnetization
      a_PoloidalMom[dit].copy(guidingCenterPoloidalMom[dit]);
      a_PoloidalMom[dit].plus(curlb_pol[dit]);
      a_PoloidalMom[dit].plus(gradM_r[dit]);
   }


}

void KineticSpecies::ParticleFlux( CFG::LevelData<CFG::FArrayBox>& a_ParticleFlux,
                                   const LevelData<FluxBox>& field  ) const

{
   m_moment_op.compute( a_ParticleFlux, *this, ParticleFluxKernel(field) );

   //Calculate flux average
   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   CFG::FluxSurface m_flux_surface(mag_geom, false);
   CFG::LevelData<CFG::FArrayBox> FluxAver_tmp;
   FluxAver_tmp.define(a_ParticleFlux);
   m_flux_surface.averageAndSpread(a_ParticleFlux, FluxAver_tmp);
   m_flux_surface.averageAndSpread(FluxAver_tmp,a_ParticleFlux);

}

void KineticSpecies::HeatFlux( CFG::LevelData<CFG::FArrayBox>& a_HeatFlux,
                                   const LevelData<FluxBox>& field,
                                   const LevelData<FArrayBox>& phi  ) const

{
   m_moment_op.compute( a_HeatFlux, *this, HeatFluxKernel(field, phi) );

   //Calculate flux average
   const CFG::MagGeom& mag_geom = m_geometry->magGeom();
   CFG::FluxSurface m_flux_surface(mag_geom, false);
   CFG::LevelData<CFG::FArrayBox> FluxAver_tmp;
   FluxAver_tmp.define(a_HeatFlux);
   m_flux_surface.averageAndSpread(a_HeatFlux, FluxAver_tmp);
   m_flux_surface.averageAndSpread(FluxAver_tmp,a_HeatFlux);

}

void KineticSpecies::parallelHeatFluxMoment( CFG::LevelData<CFG::FArrayBox>& a_parallelHeatFlux,
                                             CFG::LevelData<CFG::FArrayBox>& a_vparshift ) const
{
   m_moment_op.compute( a_parallelHeatFlux, *this, ParallelHeatFluxKernel(a_vparshift) );
}

void KineticSpecies::perpCurrentDensity( CFG::LevelData<CFG::FArrayBox>& a_PerpCurrentDensity,
                                         const LevelData<FluxBox>& field  ) const
{
   m_moment_op.compute( a_PerpCurrentDensity, *this, PerpCurrentDensityKernel(field) );
   
}

void KineticSpecies::pressureMoment( CFG::LevelData<CFG::FArrayBox>& a_pressure,
                                     CFG::LevelData<CFG::FArrayBox>& a_vparshift ) const
{
   m_moment_op.compute( a_pressure, *this, PressureKernel(a_vparshift) );
}

void KineticSpecies::energyMoment( CFG::LevelData<CFG::FArrayBox>& a_energy ) const
{
   m_moment_op.compute( a_energy, *this, EnergyKernel() );
}

void KineticSpecies::energyMoment( CFG::LevelData<CFG::FArrayBox>&  a_energy,
                                   const LevelData<FArrayBox>&      a_function ) const
{
   m_moment_op.compute( a_energy, *this, a_function, EnergyKernel() );
}

void KineticSpecies::fourthMoment( CFG::LevelData<CFG::FArrayBox>& a_fourth ) const
{
   m_moment_op.compute( a_fourth, *this, FourthMomentKernel() );
}


void
KineticSpecies::computeFSavgMaxwellian( LevelData<FArrayBox>&  a_F0 ) const
{
   
   /* Computes Vpar-unshifted FS averaged maxwellian fit. No B_star_par or J factors are included */
   const CFG::MagGeom& mag_geom( m_geometry->magGeom() );

   IntVect ghost_vect = a_F0.ghostVect();
   CFG::IntVect ghost_cfg;
   ghost_cfg[RADIAL_DIR] = ghost_vect[RADIAL_DIR];
   ghost_cfg[POLOIDAL_DIR] = ghost_vect[POLOIDAL_DIR];
   
   CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_moment_op.compute( density, *this, DensityKernel() );

   CFG::LevelData<CFG::FArrayBox> parallelMom( mag_geom.grids(), 1, ghost_cfg );
   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      parallelMom[dit].setVal(0.0);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   m_moment_op.compute( pressure, *this, PressureKernel(parallelMom) );
   

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      pressure[dit].divide(density[dit]);
   }
   
   CFG::FluxSurface fs(mag_geom);
   
   CFG::LevelData<CFG::FArrayBox> fs_average_density(mag_geom.grids(), 1, CFG::IntVect::Zero);
   fs.averageAndSpread(density, fs_average_density);
   
   CFG::LevelData<CFG::FArrayBox> fs_average_parMom(mag_geom.grids(), 1, CFG::IntVect::Zero);
   fs.averageAndSpread(parallelMom, fs_average_parMom);
   
   CFG::LevelData<CFG::FArrayBox> fs_average_temperature(mag_geom.grids(), 1, CFG::IntVect::Zero);
   fs.averageAndSpread(pressure, fs_average_temperature);
   
   
   MaxwellianKernel maxwellian(fs_average_density,
                               fs_average_temperature,
                               fs_average_parMom);

   const DisjointBoxLayout& dbl( a_F0.getBoxes() );
   LevelData<FArrayBox> dfn_no_ghost( dbl, 1, IntVect::Zero );

   maxwellian.eval(dfn_no_ghost, *this);

   for (DataIterator dit(a_F0.dataIterator()); dit.ok(); ++dit) {
     a_F0[dit].copy(dfn_no_ghost[dit]);
   }

   for (DataIterator dit(a_F0.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = m_geometry->getBlockCoordSys(dbl[dit]);
      const ProblemDomain& domain = coord_sys.domain();
      fourthOrderCellExtrapAtDomainBdry(a_F0[dit], domain, dbl[dit]);
   }
   m_geometry->fillInternalGhosts(a_F0);
   //NB: divideBstarParallel also includes fillInternalGhosts, 
   //but we included one above, just for the clarity of the code
   m_geometry->divideBStarParallel( a_F0 );

}


bool KineticSpecies::isSpecies( const string& a_name ) const
{
   if (m_name == a_name) return true;
   return false;
}


const KineticSpecies& KineticSpecies::operator=( const KineticSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      m_name = a_rhs.m_name;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_dist_func.define( a_rhs.m_dist_func );
   }
   return *this;
}


void KineticSpecies::copy( const KineticSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      m_name = a_rhs.m_name;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;

      DataIterator dit( m_dist_func.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         m_dist_func[dit].copy( a_rhs.m_dist_func[dit] );
      }
   }
}


void KineticSpecies::zeroData()
{
   DataIterator dit( m_dist_func.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      m_dist_func[dit].setVal( 0.0 );
   }
}


void KineticSpecies::addData( const KineticSpecies& a_rhs,
                              const Real a_factor )
{
   try {
      DataIterator dit( m_dist_func.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         m_dist_func[dit].plus( a_rhs.m_dist_func[dit], a_factor );
      }
   }
   catch (std::bad_cast) {
      MayDay::Error( "Invalid SpeciesModel passed to KineticSpecies::addData!" );
   }
}

bool KineticSpecies::conformsTo( const KineticSpecies& a_rhs,
                                 const bool a_include_ghost_cells ) const
{
   try {
      const LevelData<FArrayBox>& thisData = m_dist_func;
      const LevelData<FArrayBox>& rhsData = a_rhs.m_dist_func;

      const DisjointBoxLayout& thisBoxes = thisData.disjointBoxLayout();
      const DisjointBoxLayout& rhsBoxes = rhsData.disjointBoxLayout();

      bool status( true );
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( thisData.nComp() == rhsData.nComp() );

      if ( a_include_ghost_cells) {
         status &= ( thisData.ghostVect() == rhsData.ghostVect() );
      }

      return status;
   }
   catch (std::bad_cast) {
      MayDay::Error( "Invalid SpeciesModel passed to KineticSpecies::comformsTo!" );
   }
   return false;
}

RefCountedPtr<KineticSpecies>
KineticSpecies::clone( const IntVect ghostVect,
                       const bool copy_soln_data ) const
{
   RefCountedPtr<KineticSpecies> result
      = RefCountedPtr<KineticSpecies>(
              new KineticSpecies( m_name, m_mass, m_charge, m_geometry ) );

   result->m_dist_func.define( m_dist_func.disjointBoxLayout(),
                               m_dist_func.nComp(),
                               ghostVect );

   if (copy_soln_data) {
      LevelData<FArrayBox>& result_dfn = result->m_dist_func;
      DataIterator dit( result_dfn.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         result_dfn[dit].copy( m_dist_func[dit] );
      }
   }

   return result;
}


const DisjointBoxLayout& KineticSpecies::getGhostDBL() const
{
   if ( m_ghost_dbl.size() == 0 ) {

      CH_assert(m_dist_func.isDefined());
      const DisjointBoxLayout& grids = m_dist_func.disjointBoxLayout();
      const IntVect& ghosts = m_dist_func.ghostVect();

      m_ghost_dbl.deepCopy(grids);

      const PhaseCoordSys& coord_sys = m_geometry->phaseCoordSys();

      for (int block=0; block<coord_sys.numBlocks(); ++block) {
         const PhaseBlockCoordSys& block_coord_sys = (const PhaseBlockCoordSys&)(*coord_sys.getCoordSys(block));
         const ProblemDomain& domain = block_coord_sys.domain();
         const Box& domain_box = domain.domainBox();

         for (int i=0; i<m_ghost_dbl.rawPtr()->size(); ++i) {
            Box& new_box = *const_cast<Box*>(&((*m_ghost_dbl.rawPtr())[i].box));
            Box original_box = new_box;

            if ( domain_box.contains(original_box) ) {
               for (int dir=0; dir<CFG_DIM; ++dir) {
                  if ( !domain.isPeriodic(dir) ) {
                     if ( original_box.smallEnd(dir) == domain_box.smallEnd(dir) ) {
                        new_box.growLo(dir,ghosts[dir]);
                     }
                     if ( original_box.bigEnd(dir) == domain_box.bigEnd(dir) ) {
                        new_box.growHi(dir,ghosts[dir]);
                     }
                  }
               }
            }
         }
      }

      m_ghost_dbl.closeNoSort();
   }

   return m_ghost_dbl;
}


Real KineticSpecies::maxValue() const
{
   const DisjointBoxLayout& grids( m_dist_func.getBoxes() );
   Real local_maximum( -CH_BADVAL );
   DataIterator dit( m_dist_func.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      Box box( grids[dit] );
      Real box_max( m_dist_func[dit].max( box ) );
      local_maximum = Max( local_maximum, box_max );
   }

   Real maximum( local_maximum );
#ifdef CH_MPI
   MPI_Allreduce( &local_maximum, &maximum, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

   return maximum;
}

Real KineticSpecies::minValue() const
{
   const DisjointBoxLayout& grids( m_dist_func.getBoxes() );
   Real local_minimum( CH_BADVAL );
   DataIterator dit( m_dist_func.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      Box box( grids[dit] );
      Real box_min( m_dist_func[dit].min( box ) );
      if (box_min<0.0) {
         IntVect box_min_loc = m_dist_func[dit].minIndex( box );
         pout() << box_min_loc << ":  " << box_min << endl;
      }
      local_minimum = Min( local_minimum, box_min );
   }
   Real minimum(local_minimum);
#ifdef CH_MPI
   MPI_Allreduce( &local_minimum, &minimum, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
#endif
   return minimum;
}


void KineticSpecies::computeVelocity(LevelData<FluxBox>& a_velocity,
                                     const LevelData<FluxBox>& a_E_field ) const
{
   const DisjointBoxLayout& dbl( m_dist_func.getBoxes() );
   a_velocity.define( dbl, SpaceDim, IntVect::Unit );
   m_geometry->updateVelocities( a_E_field, a_velocity, PhaseGeom::FULL_VELOCITY, true );
}


void KineticSpecies::computeMappedVelocityNormals(LevelData<FluxBox>& a_velocity_normal,
                                                  const CFG::LevelData<CFG::FArrayBox>& a_E_field_cell,
                                                  const CFG::LevelData<CFG::FArrayBox>& a_phi_node,
                                                  const bool a_fourth_order_Efield,
                                                  const int  a_velocity_option) const
{
   const DisjointBoxLayout& dbl( m_dist_func.getBoxes() );
   a_velocity_normal.define( dbl, 1, IntVect::Unit );
   m_geometry->updateVelocityNormals( a_E_field_cell, a_phi_node, a_fourth_order_Efield, a_velocity_normal, a_velocity_option);
}


void KineticSpecies::computeMappedVelocity(LevelData<FluxBox>& a_velocity,
                                           const LevelData<FluxBox>& a_E_field ) const
{
   const DisjointBoxLayout& dbl( m_dist_func.getBoxes() );
   a_velocity.define( dbl, SpaceDim, IntVect::Unit );
   m_geometry->updateMappedVelocities( a_E_field, a_velocity );
}


#include "NamespaceFooter.H"
