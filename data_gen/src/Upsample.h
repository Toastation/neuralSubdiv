#include "RandomDecimation.h"

#include <vector>

#include <pmp/SurfaceMesh.h>

namespace neuralSubdiv {

    // upsample without moving vertices
	void upsample_mid_point(pmp::SurfaceMesh& mesh, pmp::SurfaceMesh& coarse);

    void ssp(int data_id, pmp::SurfaceMesh& coarse,
             pmp::SurfaceMesh& original, pmp::SurfaceMesh& mapped,
             const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos);

	void coarse_to_fine(Eigen::Array3d& bary_coord, Eigen::Array3i& bary_face, 
						const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos);

    static void inline uv_bary_coord(Eigen::Vector2d& point, Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv, Eigen::MatrixXd& bary_uv)
    {
        Eigen::Vector2d a, b, c;
        Eigen::Vector2d p0, p1, p2;
        double d00, d01, d11, d20, d21;
        double denom;
        double u, v, w;
        for (int i = 0; i < F_uv.rows(); ++i)
        {
            a = uv.row(F_uv(i, 0));
            b = uv.row(F_uv(i, 1));
            c = uv.row(F_uv(i, 2));
            p0 = b - a;
            p1 = c - a;
            p2 = point - a;
            d00 = p0.dot(p0);
            d01 = p0.dot(p1);
            d11 = p1.dot(p1);
            d20 = p2.dot(p0);
            d21 = p2.dot(p1);
            denom = d00 * d11 - d01 * d01;
            v = (d11 * d20 - d01 * d21) / denom;
            w = (d00 * d21 - d01 * d20) / denom;
            u = 1.0 - v - w;
            bary_uv(i, 0) = u;
            bary_uv(i, 1) = v;
            bary_uv(i, 2) = w;
        }
    }
}