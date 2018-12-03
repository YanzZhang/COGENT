#include "GKCollisions.H"

#include <float.h>
#include <sstream>
#include <algorithm>

#include "NamespaceHeader.H"

GKCollisions::GKCollisions( const int a_verbose )
   : m_verbose(a_verbose)
{
bool more_kinetic_species(true);
int count(0);
while (more_kinetic_species) {
   // look for all the kinetic species firstly
   std::stringstream s;
   s << "kinetic_species." << count+1;
   ParmParse ppspecies( s.str().c_str() );

   std::string species_name("Invalid");
   if (ppspecies.contains("name")) {
      ppspecies.get("name", species_name);
       m_species_name.push_back(species_name);
      typedef std::map<std::string,int>::value_type valType;
      m_species_map.insert( valType( species_name, count) );
      count++;
   }
   else {
      more_kinetic_species = false;
   }
}
//try to get all the collisional operators for all species
m_num_species=count;
count=0;
while (count<m_num_species) {
  std::stringstream ss;
   // look for all collisional operator for all the kinetic species
  ss << "kinetic_species." << count+1;
  std::string species_name("Invalid");
   ParmParse ppspecies( ss.str().c_str() );
      ppspecies.get("name", species_name);
      std::string cls_type(_CLS_NONE_);
      CLSInterface* cls(NULL);

    for(int count_bkgr(0); count_bkgr<m_num_species; count_bkgr++){
      std::stringstream cls_s;
      cls_s << "cls_" << count_bkgr+1;
      std::string cls_species(cls_s.str());
      //input sample kinetic_species.1.cls_2="Krook" means species.1 collids with species 2
      std::string species_name_bkgr(m_species_name[count_bkgr]);

      if (ppspecies.contains( cls_species )) {
         ppspecies.get( cls_species.c_str(), cls_type );
         //collisional parameters input method is CLS.species_name.cls_2


         const std::string prefix(  "CLS."+ species_name+"."+species_name_bkgr);


         if (cls_type == _CLS_KROOK_) {
            cls = new Krook( prefix, m_verbose );
         }
         else if (cls_type == _CLS_MYKROOK_) {
            cls = new MyKrook( prefix, m_verbose );
         }
         else if (cls_type == _CLS_LORENTZ_) {
            cls = new Lorentz( prefix, m_verbose );
         }
         else if (cls_type == _CLS_LINEARIZED_) {
            cls = new Linearized( prefix, m_verbose );
         }
         else if (cls_type == _CLS_FOKKERPLANCK_) {
            bool like_species(false);
            if(count==count_bkgr) like_species=true;
            cls = new FokkerPlanck( prefix, m_verbose,like_species);
         }
         else if (cls_type == _CLS_CONSDRAGDIFF_) {
            cls = new ConsDragDiff( species_name, prefix, m_verbose );
         }
         else if (cls_type == _CLS_NONE_) {
            cls = new NullCLS();
         }
         else {
            if (procID()==0) {
              cout << "Unknown collision model specified for " << species_name <<" with "<<species_name_bkgr;
              cout << ". Using none.\n";
            }
            cls_type = _CLS_NONE_;
            cls = new NullCLS();
         }
      }
      else {
         if (procID()==0) cout << "No collision model specified for " << species_name <<" with "<<species_name_bkgr<< ".\n";
         cls_type = _CLS_NONE_;
         cls = new NullCLS();
      }
      m_collision_model.push_back( cls );
      m_collision_model_name.push_back(cls_type);
      //m_collision_species_bkgr store the corresponding background species_count
      m_collision_species_bkgr.push_back(count_bkgr);
      m_collision_species.push_back(count);
      if (!procID()) {
        cout << "Collision model for " << species_name <<" with " << species_name_bkgr<< "\t"  << ":\t";
        cout << cls_type << "\n" ;
      }

    }
    count++;
  }
}


GKCollisions::~GKCollisions()
{
   for (int i(0); i<m_collision_model.size(); i++ ) {
      delete m_collision_model[i];
   }
}


CLSInterface& GKCollisions::collisionModel( const std::string& a_name,const std::string& a_name_bkgr )

{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index_a((*it).second);

   const mapIterator it_bkgr( m_species_map.find( a_name_bkgr ) );
   CH_assert(it_bkgr!=m_species_map.end());
   const int index_b((*it_bkgr).second);
   return *(m_collision_model[index_a*m_num_species+index_b]);
}


std::string GKCollisions::collisionModelName( const std::string& a_name, const std::string& a_name_bkgr )
{
  typedef std::map<std::string,int>::iterator mapIterator;
  const mapIterator it( m_species_map.find( a_name ) );
  CH_assert(it!=m_species_map.end());
  const int index_a((*it).second);

  const mapIterator it_bkgr( m_species_map.find( a_name_bkgr ) );
  CH_assert(it_bkgr!=m_species_map.end());
  const int index_b((*it_bkgr).second);
   return(m_collision_model_name[index_a*m_num_species+index_b]);
}


void GKCollisions::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const bool                   a_implicit,
                                  const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      for(int species_bkgr(0); species_bkgr<a_rhs.size(); species_bkgr++){
        KineticSpecies& rhs_species_bkgr( *(a_rhs[species_bkgr]) );
        const std::string species_name_bkgr( rhs_species_bkgr.name() );
        CLSInterface& CLS( collisionModel( species_name,species_name_bkgr ) );

        if (a_implicit) {
          CLS.evalClsRHSImplicit( a_rhs, a_soln, species,species_bkgr, a_time );
        } else {
          CLS.evalClsRHSExplicit( a_rhs, a_soln, species,species_bkgr, a_time );
        }

     }
  }
}

/*void GKCollisions::accumulateApproxRHS( KineticSpeciesPtrVect&       a_rhs,
                                        const KineticSpeciesPtrVect& a_soln,
                                        const Real                   a_time )
{
  for (int species(0); species<a_rhs.size(); species++) {
     KineticSpecies& rhs_species( *(a_rhs[species]) );
     const std::string species_name( rhs_species.name() );
     for(int species_bkgr(0); species_bkgr<a_rhs.size(); species_bkgr++){
       KineticSpecies& rhs_species_bkgr( *(a_rhs[species_bkgr]) );
       const std::string species_name_bkgr( rhs_species_bkgr.name() );

     CLSInterface& CLS( collisionModel( species_name,species_name_bkgr ) );
     CLS.evalClsApproxRHS( a_rhs, a_soln, species,species_bkgr, a_time );
   }
  }
}*/


Real GKCollisions::computeDt( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
    Real tmp = m_collision_model[species*m_num_species+species_bkgr]->computeDt(soln,species,species_bkgr);
    dt = (tmp < dt ? tmp : dt);
    count++;
     }
  }
  return (count ? dt : -1);
}

Real GKCollisions::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
  Real scale = DBL_MAX;
  int count = 0;
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
    Real tmp = m_collision_model[species*m_num_species+species_bkgr]->TimeScale(soln,species,species_bkgr);
    scale = (tmp < scale ? tmp : scale);
    count++;
     }
  }
  return (count ? scale : -1);
}

bool GKCollisions::isLinear()
{
  bool linear_flag = true;
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
      linear_flag = (linear_flag && m_collision_model[species*m_num_species+species_bkgr]->isLinear());
     }
  }
  return linear_flag;
}


int GKCollisions::precondMatrixBands()
{
  int max_bands = 0;
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
    max_bands = std::max(max_bands, m_collision_model[species*m_num_species+species_bkgr]->precondMatrixBands());
    }
  }
  return max_bands;
}


void GKCollisions::assemblePrecondMatrix( void *a_P,
                                          const KineticSpeciesPtrVect& a_soln,
                                          const GlobalDOFKineticSpeciesPtrVect& a_gdofs,
                                          const Real a_shift)
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
    KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
    const std::string         species_name_bkgr(soln_species_bkgr.name());

    CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));

    CLS.assemblePrecondMatrix(a_P,soln_species,gdofs_species,a_shift);
   }
  }
}

void GKCollisions::defineMultiPhysicsPC(std::vector<Preconditioner<GKVector,GKOps>*>& a_pc,
                                        std::vector<DOFList>&                         a_dof_list,
                                        const KineticSpeciesPtrVect&                  a_soln,
                                        const GlobalDOFKineticSpeciesPtrVect&         a_gdofs,
                                        const GKVector&                               a_x,
                                        GKOps&                                        a_ops,
                                        const std::string&                            a_out_string,
                                        const std::string&                            a_opt_string,
                                        bool                                          a_im )
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
      KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
      const std::string         species_name_bkgr(soln_species_bkgr.name());

      CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));

      CLS.defineBlockPC(a_pc, a_dof_list, a_x, a_ops, a_out_string, a_opt_string, a_im,
                      soln_species, gdofs_species, species);
     }
   }
}
void GKCollisions::updateMultiPhysicsPC(std::vector<Preconditioner<GKVector,GKOps>*>& a_pc,
                                        const KineticSpeciesPtrVect&                  a_soln,
                                        const GlobalDOFKineticSpeciesPtrVect&         a_gdofs,
                                        const Real                                    a_time,
                                        const Real                                    a_shift,
                                        const bool                                    a_im )
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
      KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
      const std::string         species_name_bkgr(soln_species_bkgr.name());

      CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));

      CLS.updateBlockPC(a_pc, soln_species, gdofs_species, a_time, a_shift, a_im, species);
     }
   }
}

void GKCollisions::preTimeStep( const KineticSpeciesPtrVect& a_soln,
                                const Real a_time,
                                const KineticSpeciesPtrVect& a_soln_physical )

{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
    KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
    const std::string         species_name_bkgr(soln_species_bkgr.name());

    CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));
    CLS.preTimeStep( a_soln, species,species_bkgr, a_time, a_soln_physical );
   }
 }
}

void GKCollisions::postTimeStage( const KineticSpeciesPtrVect& a_soln, const Real a_time, const int a_stage )
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
    KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
    const std::string         species_name_bkgr(soln_species_bkgr.name());

    CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));
    CLS.postTimeStage( a_soln, species,species_bkgr, a_time, a_stage );
   }
 }
}

#include "NamespaceFooter.H"
