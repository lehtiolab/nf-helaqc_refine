class Utils {

  public static String get_field_nr(fn, fieldname) {
    return "\$(head -n1 ${fn} | tr '\\t' '\\n' | grep -wn '^${fieldname}\$' | cut -f 1 -d':')"
  }

}
