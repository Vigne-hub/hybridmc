#ifndef INTERNALSTATE_H
#define INTERNALSTATE_H
#include <vector>
#include <fstream>

#include "vecmat3.h"
#include "molecule.h"


//#define TEST_SPHERICAL_POLAR


class DihedralState
{
    private:
        int numDihedrals;
        std::vector<double> dihedrals;
        std::vector<double> bond_angles;
        std::vector<double> bond_lengths; 
        std::vector<bool> dihedralBinary;
        

    public:
        DihedralState(): numDihedrals(0) {}
        DihedralState(int numAtoms,std::vector<Vector> &pos): numDihedrals(numAtoms-3)
        {
            //std::cout << "  Constructor: DihedralState(int numAtoms,Vector *pos): numDihedrals(numAtoms-3):" << std::endl;
            dihedrals.resize(numAtoms-3);
            bond_angles.resize(numAtoms-3);
            bond_lengths.resize(numAtoms-3);
            dihedralBinary.resize(numAtoms-3);
            formDihedrals(numAtoms, pos);
        }
        DihedralState(const Molecule &mol){
            int numAtoms = mol.getSize();
            numDihedrals = numAtoms-3;
            dihedrals.resize(numAtoms-3);
            bond_angles.resize(numAtoms-3);
            bond_lengths.resize(numAtoms-3);
            dihedralBinary.resize(numAtoms-3);
            std::vector<Vector> pos(numAtoms);
            std::vector<std::vector<double> > molPos = mol.getPositions(); // must convert to vecmat3.h format

            for (int i=0;i<numAtoms;i++) {
                pos[i].x = molPos[i][0];
                pos[i].y = molPos[i][1];
                pos[i].z = molPos[i][2];

            }
            formDihedrals(numAtoms, pos);

        }

        DihedralState(const DihedralState& other): numDihedrals(other.numDihedrals)
        {
            dihedrals.resize(other.dihedrals.size() );
            bond_angles.resize( other.bond_angles.size() );
            bond_lengths.resize(other.bond_lengths.size() );
            dihedralBinary.resize(other.dihedralBinary.size());
            for (int i=0;i<(int)other.dihedrals.size();i++)
            {
                dihedrals[i] = other.dihedrals[i];
                bond_angles[i] = other.bond_angles[i];
                bond_lengths[i] = other.bond_lengths[i];
                dihedralBinary[i] = other.dihedralBinary[i];
            }
        }
        ~DihedralState()
        {
            dihedrals.clear();
            bond_angles.clear();
            bond_lengths.clear();
            dihedralBinary.clear();
            numDihedrals = 0;
        }

        void formDihedrals(int numAtoms, std::vector<Vector> &pos)
        {
            for (int i=2;i<numAtoms-1;i++)
            {
                Vector r01 = pos[i-2] - pos[i-1];
                Vector r12 = pos[i] - pos[i-1]; // i-1 is center
                Vector r23 = pos[i+1] - pos[i]; // i is center

                double nrm01 = r01.nrm();
                double nrm12 = r12.nrm();
                double nrm23 = r23.nrm();

                Vector r01_hat = r01/nrm01;
                Vector r12_hat = r12/nrm12;
                Vector r23_hat = r23/nrm23;

                Vector z = r12_hat;
                Vector xv = z ^ r01_hat;
                Vector x = xv/xv.nrm();
                Vector y = z ^ x;

                bond_lengths[i-2] = nrm23; // length of bond from i -> i+1
                double cos_alpha = r12_hat|r23_hat; // cosine of bond angle to i+1 from i
                bond_angles[i-2] = cos_alpha;

                // Vector n1 = r01^r12; // first plane normal vector
                // Vector n2 = r23^r12; // second plane normal vector

                // Vector n1_hat = n1/n1.nrm();
                // Vector n2_hat = n2/n2.nrm();

                // Vector yt = r01_hat - (r01_hat|r12_hat) * r12_hat;
                
                double ux = r23_hat|x;
                double uy = r23_hat|y;
                double phi = 0.0;
                if (ux != 0)
                {
                    double dphi = 0.0;
                    phi = atan(uy/ux);
                    if (ux < 0.0) dphi = M_PI;

                    phi += dphi;
                }
                else{
                    phi = M_PI/2.0;
                }
                if (phi < 0.0) phi += 2.0*M_PI; // keep in range [0,2*pi]
                dihedrals[i-2] = phi/M_PI; // cos(phi) = 1 when in cis , -1 in trans

                if (fabs ( dihedrals[i-2] - 1.0) < 0.5 ) dihedralBinary[i-2] = true; // cis is true

                

#ifdef TEST_SPHERICAL_POLAR

                std::cout << "  dihedral is " << phi/M_PI << " bond angle is " << bond_angles[i-2] << std::endl;

                //  Test coordinate system for redrawing position i+1


                double part_x = r23 | x;
                double part_y = r23 | y;
                double part_z = r23 | z;

                double sin_alpha = sqrt(1.0 - cos_alpha*cos_alpha);
                double part_x_test = nrm23*sin_alpha*cos(phi);
                double part_y_test =  nrm23*sin_alpha * sin(phi);
                double part_z_test = nrm23*cos_alpha;


                Vector test_pos = part_x_test * x + part_y_test * y + part_z_test * z;

                Vector testPos = pos[i] + test_pos;
                Vector diff = pos[i+1] - testPos;

                std::cout << "  Now In formDihedrals difference vector of position " << i+1 << " is " << diff << " phi = " << phi/M_PI 
                        << " cos_alpha = " << cos_alpha << std::endl;

                double diff_x = diff|x;
                double diff_y = diff|y;
                double diff_z = diff|z;

   
                //std::cout << "     components of diff vector in local coordinates " << diff_x << " " << diff_y << " " << diff_z 
                //        << std::endl;
                std::cout << "   Decomposition of unit vector: x  " << part_x/nrm23 << "=" << part_x_test/nrm23 
                        << "  y: " << part_y/nrm23
                            << "=" << part_y_test/nrm23 << "   z: " << part_z /nrm23<< "=" 
                            << part_z_test/nrm23 << std::endl;
                
                std::cout << " ** z-hat unit vector is " << z << std::endl;
                std::cout << " ** y-hat unit vector is " << y << std::endl;
                std::cout << " ** x-hat unit vector is " << x << std::endl;
                std::cout << "  ** test Pos is " << testPos << std::endl << std::endl;

#endif


            }


        }

        std::vector<bool> &getDihedralBinary()
        {
            return dihedralBinary;
        }
        std::vector<double> getDihedralAngles ()
        {
            return dihedrals;
        }

        std::vector<double> &getBondLengths ()
        {
            return bond_lengths;
        }

        std::vector<double> &getBondAngles()
        {
            return bond_angles;
        }

        void getInternals(std::vector<double> &b, std::vector<double> &a, std::vector<double> &d)
        {
            b.resize(numDihedrals);
            a.resize(numDihedrals);
            d.resize(numDihedrals);
            for (int i=0;i<numDihedrals;i++)
            {
                b[i] = bond_lengths[i];
                a[i] = bond_angles[i];
                d[i] = dihedrals[i];
            }

        }



        std::vector<double> dihedralDifference(const DihedralState &otherSet)
        {
            std::vector<double> diff(dihedrals.size());
            if (dihedrals.size() != otherSet.dihedrals.size() )
            {
                std::cerr << "  Comparison of set of dihedrals of size " << dihedrals.size() << " with other set of size "
                        << otherSet.dihedrals.size() << " not allowed." << std::endl; 
                
            }

            
            for (int i=0;i<numDihedrals;i++) diff[i] = dihedrals[i] - otherSet.dihedrals[i];

            return diff;
        }

        void printDihedralBinary(std::ostream &out)
        {
            for (int i=0;i<numDihedrals;i++) out << dihedralBinary[i];
            out << std::endl;

        }

        void print(int scale=1)
        {
            std::cout << " Dihedral angles [";
            for (int i=0;i<numDihedrals-1;i++) std::cout << dihedrals[i]/double(scale) << ",";
            std::cout << dihedrals[numDihedrals-1]/double(scale) << "]" << std::endl;
        }

        void output(std::ofstream &out, bool frustrated = false)
        {
            for (int i=0;i<numDihedrals-1;i++) out << dihedrals[i]  << ",";
            out << dihedrals[numDihedrals-1] << "," << frustrated << std::endl;

        }

        void outputProjections(std::ofstream &out, bool frustrated = false)
        {
            for (int i=0;i<numDihedrals;i++) out << cos(M_PI*dihedrals[i])
                << "," << sin(M_PI*dihedrals[i]) << ",";
            out << frustrated << std::endl;

        }

        void outputInternalAngles(std::ofstream &out, bool frustrated = false)
        {

            for (int i=0;i<numDihedrals;i++) out << bond_angles[i] << ",";

            for (int i=0;i<numDihedrals;i++) out << cos(M_PI*dihedrals[i])
                << "," << sin(M_PI*dihedrals[i]) << ",";
            out << frustrated << std::endl;


        }

        void addAngles(const DihedralState &otherSet)
        {
             
            if (dihedrals.size() != otherSet.dihedrals.size() )
            {
                std::cerr << "  Cannot add set of dihedrals of size " << dihedrals.size() << " to other set of size "
                        << otherSet.dihedrals.size() << " not allowed." << std::endl; 
                
            }

            
            for (int i=0;i<numDihedrals;i++) 
            {
                dihedrals[i] += otherSet.dihedrals[i];
                //std::cout << " dihedral " << i << " is now " << dihedrals[i] << " added " << otherSet.dihedrals[i] << std::endl;
            }
        }
};

#endif
