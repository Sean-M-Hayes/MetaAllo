/**************************************************************************
**    Copyright 2017 Sean M. Hayes
**    This file is part of MetaAllo.
**
**    MetaAllo is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.

**    MetaAllo is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.

**    You should have received a copy of the GNU General Public License
**    along with MetaAllo.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#include "k_means.h"
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <fstream>

/////////////////////////////////////////////////////////////////////////////////

// Convenience Functions

/////////////////////////////////////////////////////////////////////////////////

void Write_Cluster_Distance(std::string Write_Path, const boost_matrix &Extrema, const std::vector<int> &Cluster_Assignment, const k_means_tree_filtering &Cluster_Fit)
{
    std::ofstream Cluster_Distance( Write_Path+"_Cluster_Distance.csv");

    Cluster_Distance<<"ID,Cluster_Assignment,";

    for(int i=0;i<Cluster_Fit.Centers.shape()[0];i++)
        if(Cluster_Fit.Centers_Map[i].size()>0)
            Cluster_Distance<<"Dist_"<<i<<",";

    Cluster_Distance<<"Period,Max_1,Max_2,Min,"<<std::endl;

    for(int i=0;i<Cluster_Assignment.size();i++)
    {
        Cluster_Distance<<i<<","<<Cluster_Assignment[i]<<",";

        for(int j=0;j<Cluster_Fit.Centers.shape()[0];j++)
            if(Cluster_Fit.Centers_Map[j].size()>0)
            {
                double Distance = 0;

                for(int k=0;k<Cluster_Fit.Centers.shape()[1];k++)
                    Distance+=pow(Extrema[i][k]-Cluster_Fit.Centers[j][k],2);

                Distance = sqrt(Distance);

                Cluster_Distance<<Distance<<",";
            }

        for(int j=0;j<Extrema.shape()[1];j++)
            Cluster_Distance<<Extrema[i][j]<<",";

        Cluster_Distance<<std::endl;
    }
}

void Write_Clustering_Data(std::string Write_Path, const boost_matrix &Extrema, const std::vector<int> &Cluster_Assignment)
{
    std::ofstream Clustering_Output( Write_Path+"_Clustering.csv");

    Clustering_Output<<"ID,Cluster,Period,Max_1,Max_2,Min,"<<std::endl;

    for(int i=0;i<Cluster_Assignment.size();i++)
    {
        Clustering_Output<<i<<","<<Cluster_Assignment[i]<<",";

        for(int j=0;j<Extrema.shape()[1];j++)
            Clustering_Output<<Extrema[i][j]<<",";

        Clustering_Output<<std::endl;
    }
}

void Write_Test_Clustering_Data(std::string Write_Path, const boost_matrix &Extrema, const std::vector<std::vector<int>> &Cluster_Assignment)
{
    std::ofstream Clustering_Output( Write_Path+"_Clustering.csv");

    Clustering_Output<<",";

    for(int i=0;i<Cluster_Assignment.size();i++)
        Clustering_Output<<"Cluster_Test_"<<i<<",";

    Clustering_Output<<"Period,Max_1,Max_2,Min"<<std::endl;

    for(int i=0;i<Extrema.shape()[0];i++)
    {
        Clustering_Output<<i<<",";

        for(int j=0;j<Cluster_Assignment.size();j++)
            Clustering_Output<<Cluster_Assignment[j][i]<<",";

        for(int j=0;j<Extrema.shape()[1];j++)
            Clustering_Output<<Extrema[i][j]<<",";

        Clustering_Output<<std::endl;
    }
}

std::vector<int> Build_Element_Sequence(int Sequence_Length, const k_means_tree_filtering &Clustering_Fit)
{
    std::vector<int> to_Return(Sequence_Length,0);

    for(int i=0;i<Clustering_Fit.Centers_Map.size();i++)
        for(int j=0;j<Clustering_Fit.Centers_Map[i].size();j++)
            to_Return[Clustering_Fit.Centers_Map[i][j]] = i;

    return to_Return;
}

/////////////////////////////////////////////////////////////////////////////////

// Normalize Data

/////////////////////////////////////////////////////////////////////////////////

Normalized_Data::Normalized_Data(const std::vector<boost_matrix> &In_Data):
    Means(In_Data[0].shape()[1],0),
    SDs(In_Data[0].shape()[1],0),
    Normalized_Means(In_Data[0].shape()[1],0),
    Normalized_SDs(In_Data[0].shape()[1],0)
{
    //Assumes all boost_matrix in In_Data are same size

    int N_Data_Points = 0;

    for(int i=0;i<In_Data.size();i++)
        N_Data_Points+=In_Data[i].shape()[0];

    Data.resize(boost::extents[N_Data_Points][In_Data[0].shape()[1]]);

    for(int i=0;i<In_Data[0].shape()[1];i++)
    {
        for(int j=0;j<In_Data.size();j++)
            for(int k=0;k<In_Data[j].shape()[0];k++)
                Means[i] += In_Data[j][k][i];

        Means[i] /= (double)N_Data_Points;
    }

    for(int i=0;i<In_Data[0].shape()[1];i++)
    {
        for(int j=0;j<In_Data.size();j++)
            for(int k=0;k<In_Data[j].shape()[0];k++)
                SDs[i] += pow(In_Data[j][k][i]-Means[i],2);

        SDs[i] /= (double)N_Data_Points;
        SDs[i] = sqrt(SDs[i]);
    }

    int Data_row = 0;

    for(int i=0;i<In_Data.size();i++)
        for(int j=0;j<In_Data[i].shape()[0];j++)
        {
            for(int k=0;k<In_Data[i].shape()[1];k++)
                if(SDs[k]>0)
                    Data[Data_row][k] = (In_Data[i][j][k]-Means[k])/SDs[k];
                else
                    Data[Data_row][k] = 0;

            Data_row++;
        }

    for(int i=0;i<SDs.size();i++)
    {
        Normalized_Means[i] = 0;
        Normalized_SDs[i] = 1;

//        if(Means[i]!=0)
//            Magnitude_CV += pow(SDs[i]/std::abs(Means[i]),2);
//        else
//            Magnitude_CV += pow(SDs[i],2);

//        Magnitude_SD += pow(SDs[i],2);
    }

//    Magnitude_CV = sqrt(Magnitude_CV);
//    Magnitude_SD = sqrt(Magnitude_SD);
}

Normalized_Data::Normalized_Data(const boost_matrix &In_Data):
    Data(boost::extents[In_Data.shape()[0]][In_Data.shape()[1]]),
    Means(In_Data.shape()[1],0),
    SDs(In_Data.shape()[1],0),
    Normalized_Means(In_Data.shape()[1],0),
    Normalized_SDs(In_Data.shape()[1],0)
{
    for(int i=0;i<In_Data.shape()[1];i++)
    {
            for(int j=0;j<In_Data.shape()[0];j++)
                Means[i] += In_Data[j][i];

        Means[i] /= (double)In_Data.shape()[0];
    }

    for(int i=0;i<In_Data.shape()[1];i++)
    {
            for(int j=0;j<In_Data.shape()[0];j++)
                SDs[i] += pow(In_Data[j][i]-Means[i],2);

        SDs[i] /= (double)In_Data.shape()[0];
        SDs[i] = sqrt(SDs[i]);
    }

    for(int i=0;i<In_Data.shape()[0];i++)
        for(int j=0;j<In_Data.shape()[1];j++)
                if(SDs[j]>0)
                    Data[i][j] = (In_Data[i][j]-Means[j])/SDs[j];
                else
                    Data[i][j] = 0;

    for(int i=0;i<SDs.size();i++)
    {
        Normalized_Means[i] = 0;
        Normalized_SDs[i] = 1;

//        if(Means[i]!=0)
//            Magnitude_CV += pow(SDs[i]/std::abs(Means[i]),2);
//        else
//            Magnitude_CV += pow(SDs[i],2);

//        Magnitude_SD += pow(SDs[i],2);
    }

//    Magnitude_CV = sqrt(Magnitude_CV);
//    Magnitude_SD = sqrt(Magnitude_SD);
}

/////////////////////////////////////////////////////////////////////////////////

// Node (for kd tree)

/////////////////////////////////////////////////////////////////////////////////

double Find_Median(std::vector<double> In_Vector)
{
    std::sort(In_Vector.begin(),In_Vector.end());

    if(In_Vector.size()%2==0)
    {
        return (In_Vector[In_Vector.size()/2-1]+In_Vector[In_Vector.size()/2])/2;
    }
    else
    {
        return In_Vector[(int)((double)In_Vector.size()/2-.5)];
    }
}
Node::Node(std::vector<int> Points, const boost_matrix &Point_List, int grow_limit, double target_cv):
    Point_Map(Points),Centroid(Point_List.shape()[1],0),Minima(Point_List.shape()[1],0),Maxima(Point_List.shape()[1],0)
{
    if(Points.size()>0)
    {
        //Calculate node stats

        for(int j=0;j<Point_List.shape()[1];j++)
        {
            Minima[j] = Point_List[Point_Map[0]][j];
            Maxima[j] = Point_List[Point_Map[0]][j];
        }

        for(int i=0;i<Points.size();i++)
        {
            for(int j=0;j<Point_List.shape()[1];j++)
            {
                Centroid[j]+=Point_List[Point_Map[i]][j];

                if(Minima[j]>Point_List[Point_Map[i]][j])
                    Minima[j]=Point_List[Point_Map[i]][j];

                if(Maxima[j]<Point_List[Point_Map[i]][j])
                    Maxima[j]=Point_List[Point_Map[i]][j];
            }

        }

        std::vector<double> Mean(Point_List.shape()[1],0);

        for(int j=0;j<Point_List.shape()[1];j++)
            Mean[j] = Centroid[j]/(double)Points.size();

        std::vector<double> SD(Point_List.shape()[1],0);

        for(int i=0;i<Points.size();i++)
            for(int j=0;j<Point_List.shape()[1];j++)
                SD[j] += pow(Point_List[Point_Map[i]][j]-Mean[j],2);

        std::vector<double> CV(Point_List.shape()[1],0);

        bool target_cv_met = true;

        for(int j=0;j<Point_List.shape()[1];j++)
        {
            SD[j] /= (double)Points.size();
            SD[j] = sqrt(SD[j]);

            CV[j] = SD[j]/Mean[j];

            if(CV[j]>target_cv)
                target_cv_met = false;
        }

//        if(target_cv_met)
//        {
//            left_child = NULL;
//            right_child = NULL;
//        }
//        else
//        {

            //Split points into two new hypershapes

            int Longest_Dimension = std::distance(CV.begin(),std::max_element(CV.begin(),CV.end()));

            std::vector<double> Splitting_Dimension;

            for(int i=0;i<Points.size();i++)
                Splitting_Dimension.push_back(Point_List[Point_Map[i]][Longest_Dimension]);

            double Split_Point = Find_Median(Splitting_Dimension);

            std::vector<int> left_points;
            std::vector<int> right_points;

            for(int i=0;i<Splitting_Dimension.size();i++)
            {
                if(Splitting_Dimension[i]>=Split_Point)
                    left_points.push_back(Point_Map[i]);
                else
                    right_points.push_back(Point_Map[i]);
            }

            //If all points have been grouped together, try splitting only those greater than the split point
            if((left_points.size()==0||right_points.size()==0)&&(left_points.size()+right_points.size())>1)
            {
                left_points.clear();
                right_points.clear();

                for(int i=0;i<Splitting_Dimension.size();i++)
                {
                    if(Splitting_Dimension[i]>Split_Point)
                        left_points.push_back(Point_Map[i]);
                    else
                        right_points.push_back(Point_Map[i]);
                }
            }

            //Grow children

            //if(left_points.size()>0&&right_points.size()>0&&!target_cv_met)
            if(left_points.size()>0&&right_points.size()>0)
            {
                if(left_points.size()>=grow_limit)
                    left_child = new Node(left_points,Point_List,grow_limit,target_cv);
                else
                    left_child = NULL;

                if(right_points.size()>=grow_limit)
                    right_child = new Node(right_points,Point_List,grow_limit,target_cv);
                else
                    right_child = NULL;
            }
            else
            {
                left_child = NULL;
                right_child = NULL;
            }
        //}
    }
    else
    {
        left_child = NULL;
        right_child = NULL;
    }
}

double Center_Distance(std::vector<double> Node_Weighted_Centroid, int Node_Size, const boost_matrix &Centers, const std::vector<std::vector<int>> Centers_Map, int Center_ID)
{
    double to_Return = 0;

    for(int i=0;i<Node_Weighted_Centroid.size();i++)
        to_Return+=pow((Node_Weighted_Centroid[i]/(double)Node_Size)-(Centers[Center_ID][i]/((double)Centers_Map[Center_ID].size()+1)),2);

    return sqrt(to_Return);
}

bool Center_is_Farther(int Test_Center, int Best_Center, std::vector<double> Node_Minima, std::vector<double> Node_Maxima, const boost_matrix &Centers, const std::vector<std::vector<int>> Centers_Map)
{
    std::vector<double> Node_Extrema;

    for(int i=0;i<Centers.shape()[1];i++)
    {
        double Direction = Centers[Test_Center][i]/((double)Centers_Map[Test_Center].size()+1)-Centers[Best_Center][i]/((double)Centers_Map[Best_Center].size()+1);

        if(Direction>0)
            Node_Extrema.push_back(Node_Maxima[i]);
        else
            Node_Extrema.push_back(Node_Minima[i]);
    }

    return Center_Distance(Node_Extrema,1,Centers,Centers_Map,Test_Center) >= Center_Distance(Node_Extrema,1,Centers,Centers_Map,Best_Center);
}

void Node::Update_Centers(std::vector<int> Candidates, boost_matrix &Centers, std::vector<std::vector<int>> &Centers_Map)
{
    std::vector<double> Candidate_Center_Distances;

    for(int i=0;i<Candidates.size();i++)
        Candidate_Center_Distances.push_back(Center_Distance(Centroid,Point_Map.size(),Centers,Centers_Map,Candidates[i]));

    int Best_Candidate = Candidates[std::distance(Candidate_Center_Distances.begin(),std::min_element(Candidate_Center_Distances.begin(),Candidate_Center_Distances.end()))];

    if(left_child==NULL&&right_child==NULL)
    {
        for(int i=0;i<Centers.shape()[1];i++)
            Centers[Best_Candidate][i]+=Centroid[i];

        for(int i=0;i<Point_Map.size();i++)
            Centers_Map[Best_Candidate].push_back(Point_Map[i]);
    }
    else
    {
        for(int i=0;i<Candidates.size();i++)
        {
            if(Candidates[i]==Best_Candidate)
                continue;

            if(Center_is_Farther(Candidates[i],Best_Candidate,Minima,Maxima,Centers,Centers_Map))
                Candidates.erase(Candidates.begin()+i);
        }
        if(Candidates.size()==1)
        {
            for(int i=0;i<Centers.shape()[1];i++)
                Centers[Best_Candidate][i]+=Centroid[i];

            for(int i=0;i<Point_Map.size();i++)
                Centers_Map[Best_Candidate].push_back(Point_Map[i]);

        }
        else
        {
            left_child->Update_Centers(Candidates,Centers,Centers_Map);
            right_child->Update_Centers(Candidates,Centers,Centers_Map);
        }
    }
}

Node::~Node()
{
    if(left_child)
        delete left_child;

    if(right_child)
        delete right_child;

}


double Inter_Cluster_Vector_Distance(std::vector<double> Vec_1, int Point, const boost_matrix &Points_List)
{
    double to_Return = 0;


    for(int i=0;i<Vec_1.size();i++)
        to_Return+=pow(Vec_1[i]-Points_List[Point][i],2);

    return sqrt(to_Return);

}

double Inter_Cluster_Vector_Distance(int Center, const boost_matrix &Centers, int Point, const boost_matrix &Points_List)
{
    double to_Return = 0;


    for(int i=0;i<Centers.shape()[1];i++)
        to_Return+=pow(Centers[Center][i]-Points_List[Point][i],2);

    return sqrt(to_Return);
}

/////////////////////////////////////////////////////////////////////////////////

// k_means_tree_filtering

/////////////////////////////////////////////////////////////////////////////////

std::vector<int> N_Unique_Samples(int min, int max, int N)
{
    if(max-min+1>=N)
    {

        //Sloppy solution, will infinte loop if max-min < N, don't use unless max-min >>> N

        boost::random::random_device true_rand;
        boost::random::mt19937 mt_generator(true_rand());
        boost::random::uniform_int_distribution< int > Select_Centers(min,max);

        std::vector<int> Initial_Points;

        for(int i=0;i<N;i++)
        {
            int Test_Center_Point = Select_Centers(mt_generator);

            while(true)
            {
                bool New_Point = true;

                for(int j=0;j<Initial_Points.size();j++)
                    if(Test_Center_Point==Initial_Points[j])
                    {
                        New_Point = false;
                        break;
                    }

                if(New_Point)
                    break;
                else
                    Test_Center_Point = Select_Centers(mt_generator);
            }

            Initial_Points.push_back(Test_Center_Point);
        }

        return Initial_Points;
    }
    else
        return std::vector<int>();
}

k_means_tree_filtering::k_means_tree_filtering(const boost_matrix &Test_Vectors, const std::vector<double> &Test_SDs, int k_Centers, double Convergence_Criteria, int Niter, Node &kd_tree_root):
    Centers(boost::extents[k_Centers][Test_Vectors.shape()[1]]),S_Dbw_Center_Fit(10),Convergence_Failed(false),Centers_Failed(true)
{
    int Testinger = Test_Vectors.shape()[0];

    if(Test_Vectors.shape()[0]>=k_Centers)
    {
        std::vector<int> Initial_Points = N_Unique_Samples(0,Test_Vectors.shape()[0]-1,k_Centers);

        std::vector<int> Candidate_Centers;

        if(Initial_Points.size()!=k_Centers)
            Centers.resize(boost::extents[Initial_Points.size()][Test_Vectors.shape()[1]]);
        else
            Centers_Failed = false;

        if(!Centers_Failed)
        {
            //Initialize Centers
            for(int i=0;i<Centers.shape()[0];i++)
            {
                Candidate_Centers.push_back(i);

                for(int j=0;j<Centers.shape()[1];j++)
                    Centers[i][j] = Test_Vectors[Initial_Points[i]][j];
            }

            //Run Centers through kd_tree until convergence found, or Niter
            for(int i=0;i<Niter;i++)
            {
                Centers_Map.clear();

                for(int i=0;i<Centers.shape()[0];i++)
                    Centers_Map.push_back(std::vector<int>());

                boost_matrix Old_Centers = Centers;

                kd_tree_root.Update_Centers(Candidate_Centers,Centers,Centers_Map);

                //update Centers from weighted, calculate distance change from previous

                double Max_Distance_Change = 0;

                for(int i=0;i<Centers.shape()[0];i++)
                {
                    double Center_Distance = 0;

                    for(int j=0;j<Centers.shape()[1];j++)
                    {
                        Centers[i][j] /= Centers_Map[i].size()+1;

                        Center_Distance+=pow(Centers[i][j]-Old_Centers[i][j],2);
                    }

                    Center_Distance = sqrt(Center_Distance);

                    if(Max_Distance_Change<Center_Distance)
                        Max_Distance_Change = Center_Distance;
                }

                if(Max_Distance_Change<=Convergence_Criteria)
                    break;

                if(i==(Niter-1))
                    Convergence_Failed = true;
            }

            //Purge any centers which have no points assigned to them
            std::vector<int> Live_Centers;
            for(int i=0;i<Centers_Map.size();i++)
                if(Centers_Map[i].size()>0)
                    Live_Centers.push_back(i);

            if(Live_Centers.size()<k_Centers)
            {
                boost_matrix Temp_Centers = Centers;
                std::vector<std::vector<int>> Temp_Centers_Map = Centers_Map;

                Centers.resize(boost::extents[Live_Centers.size()][Temp_Centers.shape()[1]]);
                Centers_Map.clear();

                for(int i=0;i<Centers.shape()[0];i++)
                {
                    for(int j=0;j<Centers.shape()[1];j++)
                        Centers[i][j] = Temp_Centers[Live_Centers[i]][j];

                    Centers_Map.push_back(Temp_Centers_Map[Live_Centers[i]]);
                }
            }

            //Check goodness of fit
            S_Dbw_Center_Fit = S_Dbw_Cluster_Goodness_of_Fit(Centers,Centers_Map,Test_Vectors,Test_SDs);
        }
    }
}

void k_means_tree_filtering::operator=(k_means_tree_filtering New_Values)
{
    S_Dbw_Center_Fit = New_Values.S_Dbw_Center_Fit;
    Centers_Map = New_Values.Centers_Map;

    Centers.resize(boost::extents[New_Values.Centers.shape()[0]][New_Values.Centers.shape()[1]]);

    for(int i=0;i<Centers.shape()[0];i++)
        for(int j=0;j<Centers.shape()[1];j++)
            Centers[i][j] = New_Values.Centers[i][j];

}

/////////////////////////////////////////////////////////////////////////////////

// Center goodness of fit measure

/////////////////////////////////////////////////////////////////////////////////

// Halkidi and Vazirgiannis, "Clustering validity assesment: Finding the optimal partitioning of a data set" in ICDM, Washington, DC, USA 2001, pp. 187-194
double S_Dbw_Cluster_Goodness_of_Fit(const boost_matrix &Centers, const std::vector<std::vector<int>> &Centers_Map, const boost_matrix &Points_List, const std::vector<double> &Point_SDs)
{
    double Total_Data_Variance = 0;

    for(int i=0;i<Point_SDs.size();i++)
    {
        //raised to the fourth because formula originally calls for sum of variance^2 (=SD^4),
        Total_Data_Variance+=pow(Point_SDs[i],4);
    }

    Total_Data_Variance = sqrt(Total_Data_Variance);

    boost_matrix Cluster_Dimension_Variance(boost::extents[Centers.shape()[0]][Centers.shape()[1]]);

    double N_Centers = Centers.shape()[0];

    for(int i=0;i<Centers.shape()[0];i++)
    {
        if(Centers_Map[i].size()==0)
        {
            N_Centers--;

            for(int j=0;j<Centers.shape()[1];j++)
                Cluster_Dimension_Variance[i][j] = 0;
        }
        else
            for(int j=0;j<Centers.shape()[1];j++)
            {
                Cluster_Dimension_Variance[i][j] = 0;

                for(int k=0;k<Centers_Map[i].size();k++)
                    Cluster_Dimension_Variance[i][j] += pow(Points_List[Centers_Map[i][k]][j] - Centers[i][j],2);

                Cluster_Dimension_Variance[i][j] /= (double) Centers_Map[i].size();
            }
    }

    if(N_Centers==0)
        return 0;

    std::vector<double> Cluster_Variance(Cluster_Dimension_Variance.shape()[0],0);

    double Total_Cluster_Variance;
    double Intra_Cluster_Density = 0;

    for(int i=0;i<Cluster_Dimension_Variance.shape()[0];i++)
    {
        std::vector<double> Cluster_Dummy;

        for(int j=0;j<Cluster_Dimension_Variance.shape()[1];j++)
        {
            Cluster_Variance[i] += pow(Cluster_Dimension_Variance[i][j],2);
            Cluster_Dummy.push_back(pow(Cluster_Dimension_Variance[i][j],2));
        }

        Cluster_Variance[i] = sqrt(Cluster_Variance[i]);

        Total_Cluster_Variance += Cluster_Variance[i];

        Intra_Cluster_Density += Cluster_Variance[i]/Total_Data_Variance;

    }

    Intra_Cluster_Density /= N_Centers;
    Total_Cluster_Variance = sqrt(Total_Cluster_Variance)/N_Centers;

    double Inter_Cluster_Density = 0;

    for(int i=0;i<Centers_Map.size();i++)
        for(int j=0;j<Centers_Map.size();j++)
        {
            if(i==j)
                continue;

            if(Centers_Map[i].size()==0||Centers_Map[j].size()==0)
                continue;

            std::vector<double> Midpoint;

            for(int k=0;k<Centers.shape()[1];k++)
                Midpoint.push_back((Centers[i][k]+Centers[j][k])/2);

            double Midpoint_Density=0;
            //double Center_i_Density=0;
            //double Center_j_Density=0;

            for(int k=0;k<Centers_Map[i].size();k++)
            {
                if(Inter_Cluster_Vector_Distance(Midpoint,Centers_Map[i][k],Points_List)<=Total_Cluster_Variance)
                    Midpoint_Density++;

                //if(Inter_Cluster_Vector_Distance(i,Centers,Centers_Map[i][k],Points_List)<=Total_Cluster_Variance)
                    //Center_i_Density++;

                //if(Inter_Cluster_Vector_Distance(j,Centers,Centers_Map[i][k],Points_List)<=Total_Cluster_Variance)
                    //Center_j_Density++;
            }

            for(int k=0;k<Centers_Map[j].size();k++)
            {
                if(Inter_Cluster_Vector_Distance(Midpoint,Centers_Map[j][k],Points_List)<=Total_Cluster_Variance)
                    Midpoint_Density++;

                //if(Inter_Cluster_Vector_Distance(i,Centers,Centers_Map[j][k],Points_List)<=Total_Cluster_Variance)
                    //Center_i_Density++;

                //if(Inter_Cluster_Vector_Distance(j,Centers,Centers_Map[j][k],Points_List)<=Total_Cluster_Variance)
                    //Center_j_Density++;
            }

            if(Midpoint_Density>0)
            {
                if(Centers_Map[i].size()>Centers_Map[j].size())
                    Inter_Cluster_Density += Midpoint_Density/((double)Centers_Map[i].size());
                else
                    Inter_Cluster_Density += Midpoint_Density/((double)Centers_Map[j].size());
            }

        }

    if(Inter_Cluster_Density>0)
        Inter_Cluster_Density /= N_Centers*(N_Centers-1);

    return Intra_Cluster_Density+Inter_Cluster_Density;
}

//////////////////////////

std::vector<std::vector<int>> Build_Fuzzy_Element_Sequence(const boost_matrix &Data, const k_means_tree_filtering &Clustering_Fit)
{
    std::vector<std::vector<int>> to_Return(Data.shape()[0],std::vector<int>());

    for(int i=0;i<Clustering_Fit.Centers_Map.size();i++)
        for(int j=0;j<Clustering_Fit.Centers_Map[i].size();j++)
        {
            to_Return[Clustering_Fit.Centers_Map[i][j]].push_back(i);

            double Min_Distance = Inter_Cluster_Vector_Distance(i,Clustering_Fit.Centers,Clustering_Fit.Centers_Map[i][j],Data);

            for(int k=0;k<Clustering_Fit.Centers.shape()[0];k++)
            {
                if(i==k)
                    continue;

                //The 2 is a little magic-numbery
                if(Inter_Cluster_Vector_Distance(k,Clustering_Fit.Centers,Clustering_Fit.Centers_Map[i][j],Data)/Min_Distance<2)
                    to_Return[Clustering_Fit.Centers_Map[i][j]].push_back(k);
            }
        }

    return to_Return;
}
