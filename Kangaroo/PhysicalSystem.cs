
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kangaroo
{
    public class PhysicalSystem
    {
        private int Iterations = 0;
        private List<Particle> m_particles = new List<Particle>();
        private double vSum = 0.0;
        public void Restart()
        {
            foreach (Particle local1 in this.m_particles)
            {
                local1.Position = local1.StartPosition;
                local1.Velocity = Vector3d.Zero;
                local1.Orientation = local1.StartOrientation;
                local1.AngularVelocity = Vector3d.Zero;
            }
            this.Iterations = 0;
        }
        public void AddParticle(Point3d p, double m)
        {
            this.m_particles.Add(new Particle(p, m));
        }
        public void ClearParticles()
        {
            this.m_particles.Clear();
            this.Iterations = 0;
        }
        public void DeleteParticle(int i)
        {
            this.m_particles.RemoveAt(i);
        }
        public int GetIterations() { return this.Iterations; }
        public int ParticleCount()
        {
            return this.m_particles.Count;
        }
        public void SetParticleList(List<Point3d> p)
        {
            this.ClearParticles();
            foreach (Point3d pointd in p)
            {
                this.m_particles.Add(new Particle(pointd, 1.0));
            }
        }
        public int FindParticleIndex(Point3d Pos, double tol, bool ByCurrent)
        {
            if (ByCurrent)
            {
                return this.m_particles.FindIndex(x => Util.OrthoClose(x.Position, Pos, tol));
            }
            return this.m_particles.FindIndex(x => Util.OrthoClose(x.StartPosition, Pos, tol));
        }
        public int FindOrientedParticleIndex(Plane P, double tol, bool ByCurrent)
        {
            return this.m_particles.FindIndex(x => (Util.OrthoClose(x.StartPosition, P.Origin, tol) && (x.StartOrientation.XAxis.IsParallelTo(P.XAxis) == 1)) && (x.StartOrientation.YAxis.IsParallelTo(P.YAxis) == 1));
        }
        public void AssignPIndex(IGoal Goal, double Tolerance)
        {
            this.AssignPIndex(Goal, Tolerance, false);
        }
        public void AssignPIndex(IGoal Goal, double Tolerance, bool ByCurrent)
        {
            Goal.PIndex = new int[Goal.PPos.Length];
            for (int i = 0; i < Goal.PPos.Length; i++)
            {
                int num2 = this.FindParticleIndex(Goal.PPos[i], Tolerance, ByCurrent);
                if (num2 == -1)
                {
                    this.AddParticle(Goal.PPos[i], 1.0);
                    Goal.PIndex[i] = this.ParticleCount() - 1;
                }
                else
                {
                    Goal.PIndex[i] = num2;
                }
                if ((Goal.InitialOrientation != null) && Goal.InitialOrientation[i].IsValid)
                {
                    Plane orientation = this.m_particles[Goal.PIndex[i]].Orientation;
                    if (!this.m_particles[Goal.PIndex[i]].Orientation.IsValid)
                    {
                        this.m_particles[Goal.PIndex[i]].Orientation = Goal.InitialOrientation[i];
                        this.m_particles[Goal.PIndex[i]].StartOrientation = Goal.InitialOrientation[i];
                    }
                    else
                    {
                        num2 = this.FindOrientedParticleIndex(Goal.InitialOrientation[i], Tolerance, ByCurrent);
                        if (num2 == -1)
                        {
                            this.AddParticle(Goal.PPos[i], 1.0);
                            Goal.PIndex[i] = this.ParticleCount() - 1;
                            this.m_particles[Goal.PIndex[i]].Orientation = Goal.InitialOrientation[i];
                            this.m_particles[Goal.PIndex[i]].StartOrientation = Goal.InitialOrientation[i];
                        }
                        else
                        {
                            Goal.PIndex[i] = num2;
                        }
                    }
                }
            }
        }
        public void SimpleStep(List<IGoal> goals)
        {
            this.SimpleStep(goals, true);
        }
        public void MomentumStep(List<IGoal> goals, double damping, int Iters)
        {
            for (int i = 0; i < Iters; i++)
            {
                foreach (Particle particle in this.m_particles)
                {
                    particle.Position += particle.Velocity;
                    if (particle.Orientation.IsValid)
                    {
                        particle.Orientation.Rotate(particle.AngularVelocity.Length, particle.AngularVelocity);
                        particle.Orientation.Origin = particle.Position;
                    }
                }
                Parallel.ForEach<IGoal>((IEnumerable<IGoal>)goals, (Action<IGoal>)(C => C.Calculate(this.m_particles)));
                foreach (IGoal goal in goals)
                {
                    for (int j = 0; j < goal.PIndex.Length; j++)
                    {
                        Particle local1 = this.m_particles[goal.PIndex[j]];
                        local1.MoveSum += (Vector3d)(goal.Move[j] * goal.Weighting[j]);
                        Particle local2 = this.m_particles[goal.PIndex[j]];
                        local2.WeightSum += goal.Weighting[j];
                        if (goal.Torque != null)
                        {
                            Particle local3 = this.m_particles[goal.PIndex[j]];
                            local3.TorqueSum += (Vector3d)(goal.Torque[j] * goal.TorqueWeighting[j]);
                            Particle local4 = this.m_particles[goal.PIndex[j]];
                            local4.TorqueWeightSum += goal.TorqueWeighting[j];
                        }
                    }
                }
                foreach (Particle particle in this.m_particles)
                {
                    if (!particle.MoveSum.IsZero)
                    {
                        Vector3d vectord = (Vector3d)(particle.MoveSum / (particle.WeightSum + particle.Mass));
                        particle.Position += vectord;
                        particle.Velocity += vectord;
                    }
                    particle.Velocity = (Vector3d)(particle.Velocity * damping);
                    if (!particle.TorqueSum.IsZero)
                    {
                        Vector3d axis = (Vector3d)(particle.TorqueSum / (particle.TorqueWeightSum + particle.Mass));
                        particle.Orientation.Rotate(axis.Length, axis);
                        particle.Orientation.Origin = particle.Position;
                        particle.AngularVelocity += axis;
                    }
                    particle.AngularVelocity = (Vector3d)(particle.AngularVelocity * damping);
                    particle.ClearForces();
                }
            }
            this.Iterations += Iters;
            this.vSum = 0.0;
            foreach (Particle particle in this.m_particles)
            {
                this.vSum += particle.Velocity.SquareLength;
            }
            this.vSum /= (double)this.ParticleCount();
        }
        public void SimpleStep(List<IGoal> goals, bool momentum)
        {
            foreach (Particle particle in this.m_particles)
            {
                if (momentum)
                {
                    particle.Position += particle.Velocity;
                }
                if (particle.Orientation.IsValid)
                {
                    if (momentum)
                    {
                        particle.Orientation.Rotate(particle.AngularVelocity.Length, particle.AngularVelocity);
                    }
                    particle.Orientation.Origin = particle.Position;
                }
            }
            Parallel.ForEach<IGoal>((IEnumerable<IGoal>)goals, (Action<IGoal>)(C => C.Calculate(this.m_particles)));
            foreach (IGoal goal in goals)
            {
                for (int i = 0; i < goal.PIndex.Length; i++)
                {
                    Particle local1 = this.m_particles[goal.PIndex[i]];
                    local1.MoveSum += (Vector3d)(goal.Move[i] * goal.Weighting[i]);
                    Particle local2 = this.m_particles[goal.PIndex[i]];
                    local2.WeightSum += goal.Weighting[i];
                    if (goal.Torque != null)
                    {
                        Particle local3 = this.m_particles[goal.PIndex[i]];
                        local3.TorqueSum += (Vector3d)(goal.Torque[i] * goal.TorqueWeighting[i]);
                        Particle local4 = this.m_particles[goal.PIndex[i]];
                        local4.TorqueWeightSum += goal.TorqueWeighting[i];
                    }
                }
            }
            foreach (Particle particle in this.m_particles)
            {
                if (!particle.MoveSum.IsZero)
                {
                    Vector3d vectord = (Vector3d)(particle.MoveSum / particle.WeightSum);
                    particle.Position += vectord;
                    if (momentum)
                    {
                        particle.Velocity += vectord;
                        if ((vectord * particle.Velocity) < 0.0)
                        {
                            particle.Velocity = (Vector3d)(particle.Velocity * 0.9);
                        }
                    }
                }
                else
                {
                    particle.Velocity = Vector3d.Zero;
                }
                if (!particle.TorqueSum.IsZero)
                {
                    Vector3d axis = (Vector3d)(particle.TorqueSum / particle.TorqueWeightSum);
                    particle.Orientation.Rotate(axis.Length, axis);
                    particle.Orientation.Origin = particle.Position;
                    if (momentum)
                    {
                        particle.AngularVelocity += axis;
                        if ((axis * particle.AngularVelocity) < 0.0)
                        {
                            particle.AngularVelocity = (Vector3d)(particle.AngularVelocity * 0.9);
                        }
                    }
                }
                else
                {
                    particle.AngularVelocity = Vector3d.Zero;
                }
                particle.ClearForces();
            }
            this.Iterations++;
            if (momentum)
            {
                this.vSum = 0.0;
                foreach (Particle particle in this.m_particles)
                {
                    this.vSum += particle.Velocity.SquareLength;
                }
                this.vSum /= (double)this.ParticleCount();
            }
        }
        public void Step(List<IGoal> goals, bool parallel, double ke)
        {
            if (parallel)
            {
                bool flag = false;
                //this.m_stopwatch.Restart();
                while (!flag)
                {
                    for (int i = 0; i < 10; i++)
                    {
                        foreach (Particle particle in this.m_particles)
                        {
                            particle.Position += particle.Velocity;
                            if (particle.Orientation.IsValid)
                            {
                                particle.Orientation.Origin = particle.Position;
                                particle.Orientation.Rotate(particle.AngularVelocity.Length, particle.AngularVelocity);
                            }
                        }
                        Parallel.ForEach<IGoal>((IEnumerable<IGoal>)goals, (Action<IGoal>)(C => C.Calculate(this.m_particles)));
                        foreach (IGoal goal in goals)
                        {
                            for (int j = 0; j < goal.PIndex.Length; j++)
                            {
                                Particle local1 = this.m_particles[goal.PIndex[j]];
                                local1.MoveSum += (Vector3d)(goal.Move[j] * goal.Weighting[j]);
                                Particle local2 = this.m_particles[goal.PIndex[j]];
                                local2.WeightSum += goal.Weighting[j];
                                if (goal.Torque != null)
                                {
                                    Particle local3 = this.m_particles[goal.PIndex[j]];
                                    local3.TorqueSum += (Vector3d)(goal.Torque[j] * goal.TorqueWeighting[j]);
                                    Particle local4 = this.m_particles[goal.PIndex[j]];
                                    local4.TorqueWeightSum += goal.TorqueWeighting[j];
                                }
                            }
                        }
                        foreach (Particle particle2 in this.m_particles)
                        {
                            if (!particle2.MoveSum.IsZero)
                            {
                                Vector3d vectord = (Vector3d)(particle2.MoveSum / particle2.WeightSum);
                                particle2.Position += vectord;
                                particle2.Velocity += vectord;
                                if ((vectord * particle2.Velocity) < 0.0)
                                {
                                    particle2.Velocity = (Vector3d)(particle2.Velocity * 0.9);
                                }
                                if (particle2.Orientation.IsValid)
                                {
                                    particle2.Orientation.Origin = particle2.Position;
                                }
                            }
                            else
                            {
                                particle2.Velocity = Vector3d.Zero;
                            }
                            if (!particle2.TorqueSum.IsZero)
                            {
                                Vector3d axis = (Vector3d)(particle2.TorqueSum / particle2.TorqueWeightSum);
                                particle2.Orientation.Origin = particle2.Position;
                                particle2.Orientation.Rotate(axis.Length, axis);
                                particle2.AngularVelocity += axis;
                                if ((axis * particle2.AngularVelocity) < 0.0)
                                {
                                    particle2.AngularVelocity = (Vector3d)(particle2.AngularVelocity * 0.9);
                                }
                            }
                            else
                            {
                                particle2.AngularVelocity = Vector3d.Zero;
                            }
                            particle2.ClearForces();
                        }
                        this.Iterations++;
                    }
                    this.vSum = 0.0;
                    foreach (Particle particle3 in this.m_particles)
                    {
                        this.vSum += particle3.Velocity.SquareLength;
                    }
                    this.vSum /= (double)this.ParticleCount();
                    /*
                    if ((this.m_stopwatch.ElapsedMilliseconds > 15L) || (this.vSum < ke))
                    {
                        flag = true;
                    }
                    */
                }
            }
            else
            {
                bool flag2 = false;
                //this.m_stopwatch.Restart();
                while (!flag2)
                {
                    for (int k = 0; k < 10; k++)
                    {
                        foreach (Particle particle4 in this.m_particles)
                        {
                            particle4.Position += particle4.Velocity;
                        }
                        using (List<IGoal>.Enumerator enumerator2 = goals.GetEnumerator())
                        {
                            while (enumerator2.MoveNext())
                            {
                                enumerator2.Current.Calculate(this.m_particles);
                            }
                        }
                        foreach (IGoal goal2 in goals)
                        {
                            for (int m = 0; m < goal2.PIndex.Length; m++)
                            {
                                Particle local5 = this.m_particles[goal2.PIndex[m]];
                                local5.MoveSum += (Vector3d)(goal2.Move[m] * goal2.Weighting[m]);
                                Particle local6 = this.m_particles[goal2.PIndex[m]];
                                local6.WeightSum += goal2.Weighting[m];
                            }
                        }
                        foreach (Particle particle5 in this.m_particles)
                        {
                            if (!particle5.MoveSum.IsZero)
                            {
                                Vector3d vectord3 = (Vector3d)(particle5.MoveSum / particle5.WeightSum);
                                particle5.Position += vectord3;
                                particle5.Velocity += vectord3;
                                if ((vectord3 * particle5.Velocity) < 0.0)
                                {
                                    particle5.Velocity = (Vector3d)(particle5.Velocity * 0.9);
                                }
                            }
                            else
                            {
                                particle5.Velocity = Vector3d.Zero;
                            }
                            particle5.ClearForces();
                        }
                        this.Iterations += 10;
                    }
                    this.vSum = 0.0;
                    foreach (Particle particle6 in this.m_particles)
                    {
                        this.vSum += particle6.Velocity.SquareLength;
                    }
                    this.vSum /= (double)this.ParticleCount();
                    /*
                    if ((this.m_stopwatch.ElapsedMilliseconds > 15L) || (this.vSum < ke))
                    {
                        flag2 = true;
                    }
                    */
                }
            }
        }
        public List<object> GetOutput(List<IGoal> goals)
        {
            List<object> list = new List<object>();
            foreach (IGoal goal in goals)
            {
                list.Add(goal.Output(this.m_particles));
            }
            return list;
        }
        /// //////////////////////////////
        List<IGoal> goals = new List<IGoal>();
        public void Start(List<Point3d> pts, List<IGoal> Goals)
        {
            ClearParticles();
            for (int i = 0; i < pts.Count; i++)
            {
                AddParticle(pts[i], 1.0);
            }
            goals.Clear();
            goals.AddRange(Goals);
            foreach (IGoal goal in goals)
            {
                if (goal.PPos != null)
                {
                    goal.PIndex = null;
                }
            }
            foreach (IGoal goal in goals)
            {
                if (goal.PIndex == null)
                {
                    this.AssignPIndex(goal, 0.001);
                }
            }
            Restart();
        }
    }
}
