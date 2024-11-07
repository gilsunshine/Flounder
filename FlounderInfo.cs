﻿using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace Flounder
{
  public class FlounderInfo : GH_AssemblyInfo
  {
    public override string Name => "Flounder Info";

    //Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => null;

    //Return a short string describing the purpose of this GHA library.
    public override string Description => "";

    public override Guid Id => new Guid("e01f8163-c05b-470e-afc8-eb980818990f");

    //Return a string identifying you or your company.
    public override string AuthorName => "";

    //Return a string representing your preferred contact details.
    public override string AuthorContact => "";

    //Return a string representing the version.  This returns the same version as the assembly.
    public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
  }
}