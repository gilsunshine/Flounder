using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;

namespace Flounder.Cutting
{
    
    public class PiecewiseDevelopable : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PiecewiseDevelopable class.
        /// </summary>
        public PiecewiseDevelopable()
          : base("PiecewiseDevelopable", "Piecewise Developable",
              "Piecewise developable approximation of triangle meshes.",
              "Flounder", "Cutting")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "mesh", "Mesh for calculation.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("int", "subdivisionSteps", "Number of subdivision steps",GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Resultant Mesh", "result", "piecewise developable mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Energy", "energy", "piecewise developable mesh", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {   
            int subdivisionSteps = 0;     

            // Check for inputs
            if(!DA.GetData(0, ref Globals.mesh)) { return; }
            if(!DA.GetData(1, ref subdivisionSteps)) { return; }

            // Get inputs
            DA.GetData(0, ref Globals.mesh);
            DA.GetData(1, ref subdivisionSteps);

            // Initialize Globals
            Globals.mesh.UnifyNormals();
            Globals.TopologyVertices = Globals.mesh.TopologyVertices;
            Globals.Faces = Globals.mesh.Faces;
            Globals.EigenVectors.Clear();
            Globals.Gradients.Clear();
            Globals.currentPoints.Clear();
            Globals.Energy = 0;

            int steps = subdivisionSteps + 1;

            for(int j = 0; j < steps; j++){

                for (int i = 0; i < Globals.TopologyVertices.Count; i++)
                {
                    Globals.EigenVectors.Add(new Vector3d(0,0,0));
                    Globals.Gradients.Add(new Vector3d(0,0,0));
                }

                foreach (var point in Globals.mesh.Vertices.ToPoint3dArray())
                {
                    Globals.currentPoints.Add(point);
                }

                // Flatten points into dense vector
                var doublePts = Vector<double>.Build.Dense(Globals.currentPoints.Count * 3);

                for (int i = 0; i <= doublePts.Count - 3; i += 3)
                {
                    doublePts[i] = Globals.currentPoints[i / 3].X;
                    doublePts[i + 1] = Globals.currentPoints[i / 3].Y;
                    doublePts[i + 2] = Globals.currentPoints[i / 3].Z;
                }
                
                // Run optimization
                OptimizeMesh(doublePts);

                //Update globals based on new mesh after optimization

                for (int i = 0; i < Globals.mesh.Vertices.Count; i++){
                    Globals.mesh.Vertices.SetVertex(i, Globals.currentPoints[i]);
                }

                if(j < steps - 1){
                    Globals.mesh = SubdivideFourToOne(Globals.mesh);
                    Globals.mesh.Weld(Math.PI);
                    Globals.mesh.UnifyNormals();
                    Globals.TopologyVertices = Globals.mesh.TopologyVertices;
                    Globals.Faces = Globals.mesh.Faces;
                    Globals.EigenVectors.Clear();
                    Globals.Gradients.Clear();
                    Globals.Energy = 0;
                    Globals.currentPoints.Clear();
                }
            }

            // Output optimized mesh and energy
            DA.SetData(0, Globals.mesh);
            DA.SetData(1, Globals.Energy);

            Globals.mesh = null;
        }

        public static class Globals
        {
            public static Mesh mesh;
            public static Rhino.Geometry.Collections.MeshTopologyVertexList TopologyVertices;
            public static List<Point3d> currentPoints = new List<Point3d>();
            public static Rhino.Geometry.Collections.MeshFaceList Faces;
            public static List<Vector3d> EigenVectors = new List<Vector3d>();
            public static List<Vector3d> Gradients= new List<Vector3d>();
            public static double Energy = 0;
        }

        public double SumEnergies(Vector<double> flatPoints){
            double sumEnergy = 0;

            object lockObject = new object();

            // Get energy for each topology vertex
            Parallel.For(0, Globals.TopologyVertices.Count, t =>{
                
                var matrix = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.Dense(3,3);

                // Faces with vertex coincident with topology vertex
                int[] faceIndices = Globals.TopologyVertices.ConnectedFaces(t);

                // Mesh vertices coincident with topology vertex
                int [] vertexIndices = Globals.TopologyVertices.MeshVertexIndices(t);
                
                for (int i = 0; i < faceIndices.Length; i++)
                {
                    
                    Point3d ptI = new Point3d();
                    Point3d ptJ = new Point3d();
                    Point3d ptK = new Point3d();

                    // Get face vertices
                    Point3d A = Globals.currentPoints[Globals.Faces[faceIndices[i]].A];
                    Point3d B = Globals.currentPoints[Globals.Faces[faceIndices[i]].B];
                    Point3d C = Globals.currentPoints[Globals.Faces[faceIndices[i]].C];

                    // Order face vertices such that current topology vertex is coindent with ptI
                    if(vertexIndices.Contains(Globals.Faces[faceIndices[i]].A)){ 
                        ptI = A;
                        ptJ = B;
                        ptK = C;
                    }
                    else if(vertexIndices.Contains(Globals.Faces[faceIndices[i]].B)){   
                        ptI = B;
                        ptJ = C;
                        ptK = A;
                    }
                    else if(vertexIndices.Contains(Globals.Faces[faceIndices[i]].C)){  
                        ptI = C;
                        ptJ = A;
                        ptK = B;
                    }

                    // Vectors pointing away from ptI on face edges ij and ik
                    Vector3d vij = ptJ - ptI;
                    Vector3d vik = ptK - ptI;
                    
                    vij.Unitize();
                    vik.Unitize();

                    // Calculate face normal
                    Vector3d normal = Vector3d.CrossProduct(vij, vik);
                    // I don't think this is actually necessary since vij and vik are unitized
                    normal.Unitize();

                    var normalVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);
                    normalVector[0] = normal.X;
                    normalVector[1] = normal.Y;
                    normalVector[2] = normal.Z;
                    
                    // Angle between ij and ik
                    double theta = Vector3d.VectorAngle(vij, vik);
                    
                    // Vector representing a plane perpendicular to the normal (coplanar to the face)
                    var Tmatrix = MathNet.Numerics.LinearAlgebra.Vector<double>.OuterProduct(normalVector, normalVector);
                    // Scale vector based on angle around ptI
                    Tmatrix = Tmatrix.Multiply(theta);
                    // Sum these planes around each topology vertex
                    matrix += Tmatrix;
                    
                }

                // Get the smallest eigen vector of the matrix
                var decomp = matrix.Evd();
                var eigenVectors = decomp.EigenVectors;
                var minEigenvector = eigenVectors.Column(0);
                var normalizedEigenvector = minEigenvector.Normalize(2);
                Vector3d eigenVector = new Vector3d(normalizedEigenvector[0], normalizedEigenvector[1], normalizedEigenvector[2]);
                Globals.EigenVectors[t] = eigenVector;
                double e = decomp.D[0,0];

                lock(lockObject){
                    sumEnergy += e;
                }
            });

            Globals.Energy = sumEnergy;
            return sumEnergy;
        }

        public Vector<double> AllGradients(Vector<double> flatPoints){

            for (int i = 0; i < Globals.currentPoints.Count; i++){
                Globals.currentPoints[i] = new Point3d(flatPoints[i * 3], flatPoints[i * 3 + 1], flatPoints[i * 3 + 2]);
            }
            
            // Calculate gradient for each topology vertex
            Parallel.For(0, Globals.TopologyVertices.Count, n =>
                {
                var matrix = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.Dense(3,3);
                int[] faceIndices = Globals.TopologyVertices.ConnectedFaces(n);
                int [] vertexIndices = Globals.TopologyVertices.MeshVertexIndices(n);

                Globals.Gradients[n] = new Vector3d(0,0,0);

                for (int i = 0; i < faceIndices.Length; i++)
                {
                
                    Point3d ptI = new Point3d();
                    Point3d ptJ = new Point3d();
                    Point3d ptK = new Point3d();

                    // Get face vertices
                    Point3d A = Globals.currentPoints[Globals.Faces[faceIndices[i]].A];
                    Point3d B = Globals.currentPoints[Globals.Faces[faceIndices[i]].B];
                    Point3d C = Globals.currentPoints[Globals.Faces[faceIndices[i]].C];

                    
                    int[] ptOrder = new int[3];
                    
                    // Order face vertices such that current topology vertex is coindent with ptI
                    if(vertexIndices.Contains(Globals.Faces[faceIndices[i]].A)){ 

                        ptI = A;
                        ptJ = B;
                        ptK = C;

                        ptOrder[0] = Globals.Faces[faceIndices[i]].A;
                        ptOrder[1] = Globals.Faces[faceIndices[i]].B;
                        ptOrder[2] = Globals.Faces[faceIndices[i]].C;

                    }
                    else if(vertexIndices.Contains(Globals.Faces[faceIndices[i]].B)){   

                        ptI = B;
                        ptJ = C;
                        ptK = A;

                        ptOrder[0] = Globals.Faces[faceIndices[i]].B;
                        ptOrder[1] = Globals.Faces[faceIndices[i]].C;
                        ptOrder[2] = Globals.Faces[faceIndices[i]].A;

                    }
                    else if(vertexIndices.Contains(Globals.Faces[faceIndices[i]].C)){  

                        ptI = C;
                        ptJ = A;
                        ptK = B;

                        ptOrder[0] = Globals.Faces[faceIndices[i]].C;
                        ptOrder[1] = Globals.Faces[faceIndices[i]].A;
                        ptOrder[2] = Globals.Faces[faceIndices[i]].B;

                    }
                
                    Vector3d vjk = ptK - ptJ;
                    Vector3d vij =  ptJ - ptI;
                    Vector3d vki = ptI - ptK;
                    Vector3d vji =  ptI - ptJ;
                    Vector3d vik = ptK - ptI;

                    double epsilon = 1e-10;

                    Vector3d normal = Vector3d.CrossProduct(vij, vik);

                    normal.Unitize();

                    // Calculate area safely
                    double area = Vector3d.CrossProduct(vij, vik).Length / 2.0;
                    double safeArea = Math.Max(area, epsilon);

                    var normalVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);
                    normalVector[0] = normal.X;
                    normalVector[1] = normal.Y;
                    normalVector[2] = normal.Z;

                    Vector3d cross = Vector3d.CrossProduct(normal, vjk);

                    var crossVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);
                    crossVector[0] = cross.X;
                    crossVector[1] = cross.Y;
                    crossVector[2] = cross.Z;

                    var dNi = MathNet.Numerics.LinearAlgebra.Vector<double>.OuterProduct(normalVector, crossVector) / safeArea;
                    
                    cross = Vector3d.CrossProduct(normal, vki);
                    crossVector[0] = cross.X;
                    crossVector[1] = cross.Y;
                    crossVector[2] = cross.Z;

                    var dNj = MathNet.Numerics.LinearAlgebra.Vector<double>.OuterProduct(normalVector, crossVector) / safeArea;

                    cross = Vector3d.CrossProduct(normal, vij);
                    crossVector[0] = cross.X;
                    crossVector[1] = cross.Y;
                    crossVector[2] = cross.Z;
                    
                    var dNk = MathNet.Numerics.LinearAlgebra.Vector<double>.OuterProduct(normalVector, crossVector) / safeArea;

                    vjk.Unitize();
                    vij.Unitize();
                    vki.Unitize();
                    vji.Unitize();
                    vik.Unitize();

                    Vector3d dTj = Vector3d.CrossProduct(normal, vji);
                    Vector3d dTk = Vector3d.CrossProduct(normal, vik);
                    Vector3d dTi = -dTj - dTk;

                    double theta = Vector3d.VectorAngle(vij, vik);

                    var MNeigenVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);
                    MNeigenVector[0] = Globals.EigenVectors[n][0];
                    MNeigenVector[1] = Globals.EigenVectors[n][1];
                    MNeigenVector[2] = Globals.EigenVectors[n][2];

                    var dNiVector = dNi * MNeigenVector;
                    var dNjVector = dNj * MNeigenVector;
                    var dNkVector = dNk * MNeigenVector;

                    var MNdTi = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);
                    var MNdTj = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);
                    var MNdTk = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(3);

                    for(int t = 0; t < 3; t++){
                        MNdTi[t] = dTi[t];
                        MNdTj[t] = dTj[t];
                        MNdTk[t] = dTk[t];
                    }

                    double xDotN = Globals.EigenVectors[n] * normal;

                    // Lambdas calculation
                    var dLi = MNdTi * Math.Pow(xDotN, 2.0) + dNiVector * 2 * theta * xDotN;
                    var dLj = MNdTj * Math.Pow(xDotN, 2.0) + dNjVector * 2 * theta * xDotN;
                    var dLk = MNdTk * Math.Pow(xDotN, 2.0) + dNkVector * 2 * theta * xDotN;

                    var v3DdLi = new Vector3d();
                    var v3DdLj = new Vector3d();
                    var v3DdLk = new Vector3d();

                    for(int t = 0; t < 3; t++){
                        v3DdLi[t] = dLi[t];
                        v3DdLj[t] = dLj[t];
                        v3DdLk[t] = dLk[t];
                    }

                    Globals.Gradients[n] += v3DdLi;
                    int jTop = Globals.TopologyVertices.TopologyVertexIndex(ptOrder[1]);
                    Globals.Gradients[jTop] += v3DdLj;
                    int kTop = Globals.TopologyVertices.TopologyVertexIndex(ptOrder[2]);
                    Globals.Gradients[kTop] += v3DdLk;
                }
            });
            
            Vector<double> allGradients = Vector<double>.Build.Dense(Globals.Gradients.Count * 3);

            for (int i = 0; i < Globals.Gradients.Count; i++)
            {
                allGradients[3*i] = -Globals.Gradients[i].X;
                allGradients[3*i + 1] = -Globals.Gradients[i].Y;
                allGradients[3*i + 2] = -Globals.Gradients[i].Z;

            }
            return allGradients;
                
        }

        Mesh SubdivideFourToOne (Mesh mesh){

            // Retrieve input mesh vertices and faces
            Rhino.Geometry.Collections.MeshVertexList vertices = mesh.Vertices;
            Rhino.Geometry.Collections.MeshFaceList faces = mesh.Faces;
            
            // Create new mesh to hold subdivided input mesh
            Mesh subdivisionMesh = new Mesh();

            for (int i = 0; i < faces.Count; i++) {
            
                MeshFace face = faces[i];

                // Get face vertices from input mesh
                Point3d vertexA = vertices[face.A];
                Point3d vertexB = vertices[face.B];
                Point3d vertexC = vertices[face.C];

                // Get edge midpoints for each face
                Point3d midpoint1 = Midpoint(vertexA, vertexB);
                Point3d midpoint2 = Midpoint(vertexB, vertexC);
                Point3d midpoint3 = Midpoint(vertexC, vertexA);

                // Add original vertices and midpoint vertices to new mesh
                subdivisionMesh.Vertices.Add(vertexA); //0
                subdivisionMesh.Vertices.Add(vertexB); //1
                subdivisionMesh.Vertices.Add(vertexC); //2
                subdivisionMesh.Vertices.Add(midpoint1); //3
                subdivisionMesh.Vertices.Add(midpoint2); //4
                subdivisionMesh.Vertices.Add(midpoint3); //5

                int increment = i * 6;

                // Add subdivided faces to new mesh
                subdivisionMesh.Faces.AddFace(increment + 0, increment + 3, increment + 5);
                subdivisionMesh.Faces.AddFace(increment + 3, increment + 1, increment + 4);
                subdivisionMesh.Faces.AddFace(increment + 4, increment + 2, increment + 5);
                subdivisionMesh.Faces.AddFace(increment + 3, increment + 4, increment + 5);

            }

            subdivisionMesh.Normals.ComputeNormals();
            subdivisionMesh.Compact();

            return subdivisionMesh;

        }

        Point3d Midpoint (Point3d pointA, Point3d pointB){

            // Calculate midpoint x, y and z values
            double x = (pointA.X + pointB.X) / 2;
            double y = (pointA.Y + pointB.Y) / 2;
            double z = (pointA.Z + pointB.Z) / 2;

            Point3d midpoint = new Point3d(x, y, z);

            return midpoint;

        }

        public void OptimizeMesh(Vector<double> initialVertices, int restartCount = 0)
        {
            // Termination criteria to prevent infinite loop
            int maxRestarts = 6;
            double energyChangeThreshold = 1e-5;

            if (restartCount > maxRestarts)
            {
                Console.WriteLine("Maximum restart limit reached. Terminating optimization.");
                return;
            }

            // Define the objective function
            Func<Vector<double>, double> objectiveFunction = SumEnergies;

            // Define the gradient function
            Func<Vector<double>, Vector<double>> gradientFunction = AllGradients;

            // Create the L-BFGS Minimizer (MathNet Numerics)
            var minimizer = new LimitedMemoryBfgsMinimizer(
                gradientTolerance: 1e-6,
                parameterTolerance: 1e-6,
                functionProgressTolerance: 1e-6,
                memory: 20,
                maximumIterations: 200);

            // Wrap the functions in an ObjectiveFunction
            var objective = ObjectiveFunction.Gradient(objectiveFunction, gradientFunction);

            try
            {
                double initialEnergy = objectiveFunction(initialVertices);

                // Perform the minimization
                var result = minimizer.FindMinimum(objective, initialVertices);

                double finalEnergy = result.FunctionInfoAtMinimum.Value;

                Console.WriteLine($"Optimization completed. Reason for exit: {result.ReasonForExit}");
                Console.WriteLine($"Initial energy: {initialEnergy}, Final energy: {finalEnergy}");

                // Check if energy reduction is large enough
                if (result.ReasonForExit == ExitCondition.AbsoluteGradient && Math.Abs(initialEnergy - finalEnergy) > energyChangeThreshold)
                {
                    Console.WriteLine($"Restarting optimization to further minimize energy. Restart count: {restartCount + 1}");
                    OptimizeMesh(result.MinimizingPoint, restartCount + 1);
                }

                // Update the mesh with the final vertex positions
                var optimizedVertices = result.MinimizingPoint;
                for (int i = 0; i < Globals.mesh.Vertices.Count; i++)
                {
                    Globals.mesh.Vertices.SetVertex(i,
                        (float)optimizedVertices[i * 3],
                        (float)optimizedVertices[i * 3 + 1],
                        (float)optimizedVertices[i * 3 + 2]);
                }
            }
            catch (InvalidOperationException ex)
            {
                // Handle optimization issues
                Console.WriteLine($"Optimization error: {ex.Message}");
                Console.WriteLine($"Restarting optimization due to exception. Restart count: {restartCount + 1}");
                if (restartCount < maxRestarts)
                {
                    // Perturb the vertices to try to get out of problematic configuration
                    Random rand = new Random();
                    for (int i = 0; i < initialVertices.Count; i++)
                    {
                        initialVertices[i] += (rand.NextDouble() - 0.5) * 1e-4;
                    }

                    OptimizeMesh(initialVertices, restartCount + 1);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Unhandled exception: {ex.Message}");
            }
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("80D0404D-1B13-49CB-99FB-986C9A5D509D"); }
        }
    }
}